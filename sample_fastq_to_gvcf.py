#!/usr/bin/env python3

import argparse
import boto3
import botocore
import os.path
import logging
import subprocess
import os

align_cmd = '''
bwa mem -t {threads} -R '{read_group}' {ref} {fq1} {fq2} | 
samblaster | 
samtools view -b -u /dev/stdin |
samtools sort -@ {threads} -m {mem} -O BAM -o {out} /dev/stdin
'''

sentieon_align_cmd = '''
bwa mem -t {threads} -R '{read_group}' {ref} {fq1} {fq2} |
sentieon util sort --sam2bam -o {out} -i -
'''

sentieon_dedup_cmd = '''
sentieon driver -t {threads} -i {bam} --algo LocusCollector \
--fun score_info {score}
sentieon driver -t {threads} -i {bam} --algo Dedup --rmdup \
--score_info {score} --metrics {metrics} {out}
'''

index_cmd = '''
samtools index {bam}
'''

call_vars_cmd = '''
java -Xmx{mem} -jar {gatk} -T HaplotypeCaller -R {ref} \
{input} -o {out} --emitRefConfidence GVCF -L {chrom} \
--variant_index_type LINEAR --variant_index_parameter 128000 \
-G StandardAnnotation -A AlleleBalance -A TandemRepeatAnnotator \
-A ClippingRankSumTest -A GCContent -A MappingQualityZero \
-A SpanningDeletions -A StrandOddsRatio -A AlleleBalanceBySample 
'''

sentieon_call_cmd = '''
sentieon driver -t {threads} -r {ref} {input} --algo Haplotyper \
--emit_mode GVCF --annotation AlleleBalance \
--annotation TandemRepeatAnnotator --annotation ClippingRankSumTest \
--annotation GCContent --annotation MappingQualityZero \
--annotation SpanningDeletions --annotation StrandOddsRatio \
--annotation AlleleBalanceBySample {out}
'''

gvcf_concat_cmd = '''
java -Xmx{mem} -cp {gatk} org.broadinstitute.gatk.tools.CatVariants \
-R {ref} -out {out} {input} --assumeSorted \
--variant_index_type LINEAR --variant_index_parameter 128000
'''

chroms = [str(x) for x in range(1, 23)] + ['X', 'Y', 'MT']

def download_and_align(s3, bucket, fq1, fq2, sample, read_group_id, threads, ref, mem=None, sentieon=False):
    fq1_local = "/ephemeral/" + os.path.basename(fq1)
    fq2_local = "/ephemeral/" + os.path.basename(fq2)
    
    logging.info("Downloading {} to {}".format(bucket + '/' + fq1, fq1_local))
    s3.Object(bucket, fq1).download_file(fq1_local)
    logging.info("Downloading {} to {}".format(bucket + '/' + fq2, fq2_local))
    s3.Object(bucket, fq2).download_file(fq2_local)

    read_group = r"@RG\tID:{}\tSM:{}".format(read_group_id, sample)

    if sentieon:
        out = "/ephemeral/{read_group}_dup.bam".format(read_group=read_group_id)
        cmd = sentieon_align_cmd
        cmd = cmd.format(
            threads = threads,
            read_group = read_group,
            ref = ref,
            fq1 = fq1_local,
            fq2 = fq2_local,
            out = out)
    else:
        out = "/ephemeral/{read_group}_sorted.bam".format(read_group=read_group_id)
        cmd = align_cmd
        cmd = cmd.format(
            threads = threads,
            read_group = read_group,
            ref = ref,
            fq1 = fq1_local,
            fq2 = fq2_local,
            mem = mem,
            out = out)

    logging.info("Running alignment: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    logging.info("Removing {} and {}".format(fq1_local, fq2_local))
    os.remove(fq1_local)
    os.remove(fq2_local)

    return out

def dedup_bam(bam, threads, read_group_id):
    out = "/ephemeral/{read_group}_sorted.bam".format(read_group=read_group_id)

    cmd = sentieon_dedup_cmd
    cmd = cmd.format(
        threads = threads,
        bam = bam,
        score = "/ephemeral/{read_group}_score.txt".format(read_group=read_group_id),
        metrics = "/ephemeral/{read_group}_metrics.txt".format(read_group=read_group_id),
        out = out)

    logging.info("Running duplicate removal: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    logging.info("Removing {}".format(bam))
    os.remove(bam)

    return out

def call_vars(bams, sample_name, ref, chrom=None, threads=None, mem=None, gatk=None, sentieon=False):
    if sentieon:
        input_files = " -i " + " -i ".join(bams)
    else:
        input_files = " -I " + " -I ".join(bams)

    if chrom:
        output_vars = "/ephemeral/" + sample_name + '_' + chrom + ".g.vcf.gz"
    else:
        output_vars = "/ephemeral/" + sample_name + ".g.vcf.gz"

    if sentieon:
        cmd = sentieon_call_cmd
        cmd = cmd.format(
            threads = threads,
            ref = ref,
            input = input_files,
            out = output_vars)
    else:
        cmd = call_vars_cmd
        cmd = cmd.format(
            mem = mem,
            gatk = gatk,
            ref = ref,
            input = input_files,
            out = output_vars,
            chrom = chrom)

    logging.info("Running variant calling: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    return output_vars

def index_bam(bam, threads):
    cmd = index_cmd
    cmd = cmd.format(threads=threads, bam=bam)

    logging.info("Indexing {}".format(bam))
    subprocess.check_call(cmd, shell=True)
    return

def concat_gvcfs(gvcfs, sample_name, ref, mem, gatk):
    cat_input = " -V " + " -V ".join(gvcfs)
    cat_output = "/ephemeral/" + sample_name + ".g.vcf.gz"

    cmd = gvcf_concat_cmd
    cmd = cmd.format(
        mem = mem,
        gatk = gatk,
        ref = ref,
        input = cat_input,
        out = cat_output)

    logging.info("Running cat variants: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    return cat_output

def process_args():
    parser = argparse.ArgumentParser(description="Process fastq files from a single sample to gvcf files. Initial files and final results are uploaded to S3")
    parser.add_argument("--threads", type=int, default=8, help="The number of alignment and indexing threads")
    parser.add_argument("--sort_mem", default="128M", help="Memory to use when sorting the alignment")
    parser.add_argument("--call_vars_mem", default="3g", help="Memory to use when calling variants")
    parser.add_argument("--gatk", default="/usr/local/bin/GenomeAnalysisTK.jar", help="The GATK .jar file")
    parser.add_argument("--bam_key", default="1000genomes/BAM/{sample}/{run}.bam", help="The S3 destination for temporary BAM files")
    parser.add_argument("--gvcf_key", default="1000genomes/gVCF/{sample}/{sample}_{chrom}.g.vcf.gz", help="The S3 destination for temporary gVCF files")
    parser.add_argument("--sentieon", action="store_true", help="Run the analysis with Sentieon's tools")
    parser.add_argument("reference", help="The reference genome")
    parser.add_argument("access_key", help="AWS access key")
    parser.add_argument("secret_key", help="AWS secret key")
    parser.add_argument("upload_location", help="The S3 location for file upload")
    parser.add_argument("sample_name", help="The name of the sample")
    parser.add_argument("input_fastq", nargs='+', help="The input fastq files in s3")
    return parser.parse_args()

def main(args):
    log_format = "%(asctime)s %(filename)s::%(funcName)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO,
                        format=log_format)

    logging.info("Starting analysis")

    if not args:
        args = process_args()

    session = boto3.session.Session(aws_access_key_id=args.access_key, aws_secret_access_key=args.secret_key)
    s3 = session.resource("s3")
    
    bucket_end = args.upload_location.find('/')
    bucket = args.upload_location[:bucket_end]
    if bucket.startswith("s3://"):
        bucket = bucket[5:]
    key = args.upload_location[bucket_end + 1:]

    bams = []
    s3_bams = []
    for fastq in args.input_fastq:
        fq_bucket_end = fastq.find('/')
        fq_bucket = fastq[:fq_bucket_end]
        if fq_bucket.startswith("s3://"):
            fq_bucket = fq_bucket[5:]
        fq1 = fastq[fq_bucket_end + 1:]
        fq2 = fq1[:-15] + "2.filt.fastq.gz"
        read_group_id = os.path.basename(fq1)[:-16]

        # Get alignments #
        bam_key = args.bam_key.format(sample=args.sample_name, run=read_group_id)
        s3_bams.append(bam_key)
        bam_in_s3 = False
        s3_bam = s3.Object(bucket, bam_key)
        try:
            s3_bam.load()
        except botocore.exceptions.ClientError as e:
            if e.response["Error"]["Code"] == "404":
                pass
            else:
                raise e
        else:
            bam_in_s3 = True

        if bam_in_s3:
            # Alignments are already present in s3, just download #
            bam_out = "/ephemeral/{read_group}_sorted.bam".format(read_group=read_group_id)
            logging.info("Downloading {} to {}".format(bucket + '/' + bam_key, bam_out))
            s3_bam.download_file(bam_out)
            logging.info("Downloading {} to {}".format(bucket + '/' + bam_key + ".bai", bam_out + ".bai"))
            s3.Object(bucket, bam_key + ".bai").download_file(bam_out + ".bai")
            bams.append(bam_out)
        else:
            # Make the alignments from the fastq #
            next_bam = download_and_align(s3, fq_bucket, fq1, fq2, args.sample_name, read_group_id, args.threads, args.reference, args.sort_mem, args.sentieon)
            if args.sentieon:
                previous_bam = next_bam
                next_bam = dedup_bam(previous_bam, args.threads, read_group_id)
                logging.info("Removing bam {} and index {}".format(previous_bam, previous_bam + ".bai"))
                try:
                    os.remove(previous_bam)
                except FileNotFoundError:
                    pass
                try:
                    os.remove(previous_bam + ".bai")
                except FileNotFoundError:
                    pass
            else:
                index_bam(next_bam, args.threads)
            bams.append(next_bam)
            # Upload intermediate files to s3 #
            logging.info("Uploading {} to {}".format(next_bam, bucket + '/' + bam_key))
            s3.meta.client.upload_file(next_bam, bucket, bam_key)
            logging.info("Uploading {} to {}".format(next_bam + ".bai", bucket + '/' + bam_key + ".bai"))
            s3.meta.client.upload_file(next_bam + ".bai", bucket, bam_key + ".bai")

    # Call variants on the bam files #
    gvcfs = []
    s3_gvcfs = []
    if args.sentieon:
        concat_gvcf = call_vars(bams, args.sample_name, args.reference, threads=args.threads, sentieon=True)
        for bam in bams:
            logging.info("Removing bam {} and index {}".format(bam, bam + ".bai"))
            os.remove(bam)
            os.remove(bam + ".bai")
        for tmp_file in os.listdir("/ephemeral/"):
            if tmp_file.endswith("_score.txt") or tmp_file.endswith("_metrics.txt"):
                tmp_file = "/ephemeral/" + tmp_file
                logging.info("Removing temporary file {}".format(tmp_file))
                os.remove(tmp_file)

    else:
        for chrom in chroms:
            # Check if the gvcf subset is in s3 #
            gvcf_key = args.gvcf_key.format(sample=args.sample_name, chrom=chrom)
            s3_gvcfs.append(gvcf_key)
            gvcf_in_s3 = False
            s3_gvcf = s3.Object(bucket, gvcf_key)
            try:
                s3_gvcf.load()
            except botocore.exceptions.ClientError as e:
                if e.response["Error"]["Code"] == "404":
                    pass
                else:
                    raise e
            else:
                gvcf_in_s3 = True

            if gvcf_in_s3:
                # Download the subsets #
                gvcf_out = "/ephemeral/" + args.sample_name + '_' + chrom + ".g.vcf.gz"
                logging.info("Downloading {} to {}".format(bucket + '/' + gvcf_key, gvcf_out))
                s3_gvcf.download_file(gvcf_out)
                logging.info("Downloading {} to {}".format(bucket + '/' + gvcf_key + ".tbi", gvcf_out + ".tbi"))
                s3.Object(bucket, gvcf_key + ".tbi").download_file(gvcf_out + ".tbi")
                gvcfs.append(gvcf_out)
            else:        
                next_gvcf = call_vars(bams, args.sample_name, args.reference, chrom=chrom, mem=args.call_vars_mem, gatk=args.gatk)
                gvcf_idx = next_gvcf + ".tbi"
                gvcfs.append(next_gvcf)

                # Upload the gvcf subset to s3 #
                logging.info("Upload {} to {}".format(next_gvcf, bucket + '/' + gvcf_key))
                s3.meta.client.upload_file(next_gvcf, bucket, gvcf_key)
                logging.info("Upload {} to {}".format(gvcf_idx, bucket + '/' + gvcf_key + ".tbi"))
                s3.meta.client.upload_file(gvcf_idx, bucket, gvcf_key + ".tbi")

        # Remove the BAM files #
        for bam in bams:
            logging.info("Removing bam {} and index {}".format(bam, bam + ".bai"))
            os.remove(bam)
            os.remove(bam + ".bai")

        # Concatenate the gVCF subsets #
        concat_gvcf = concat_gvcfs(gvcfs, args.sample_name, args.reference, args.call_vars_mem, args.gatk)
        for gvcf in gvcfs:
            logging.info("Removing gvcf and index {}".format(gvcf, gvcf + ".tbi"))
            os.remove(gvcf)
            os.remove(gvcf + ".tbi")

    # Upload the GVCF file #
    logging.info("Uploading {} to {}".format(concat_gvcf, bucket + '/' + key))
    s3.meta.client.upload_file(concat_gvcf, bucket, key)
    logging.info("Uploading {} to {}".format(concat_gvcf + ".tbi", bucket + '/' + key + ".tbi"))
    s3.meta.client.upload_file(concat_gvcf + ".tbi", bucket, key + ".tbi")

    logging.info("Removing gvcf file: {}".format(concat_gvcf))
    os.remove(concat_gvcf)
    logging.info("Removing index file: {}".format(concat_gvcf + ".tbi"))
    os.remove(concat_gvcf + ".tbi")

    # Clean up the intermediate files in the s3 bucket #
    logging.info("Cleaning up the s3 bucket")
    to_remove = s3_bams + [x + ".bai" for x in s3_bams] + s3_gvcfs + [x + ".tbi" for x in s3_gvcfs]
    response = s3.meta.client.delete_objects(
        Bucket = bucket,
        Delete = { 'Objects': [ {'Key':x} for x in to_remove] }
    )
    if "Errors" in response:
        logging.warning("The deletion operation returned errors: {}".format(str(response)))

    logging.info("Analsis finished")

if __name__ == "__main__":
    main(None)
