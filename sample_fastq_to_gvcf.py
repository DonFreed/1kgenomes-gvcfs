#!/usr/bin/env python3

import argparse
import boto3
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

index_cmd = '''
samtools index {bam}
'''

call_vars_cmd = '''
java -Xmx{mem} -jar {gatk} -T HaplotypeCaller -R {ref} \
{input} -o {out} --emitRefConfidence GVCF \
--variant_index_type LINEAR --variant_index_parameter 128000 \
-G StandardAnnotation -A AlleleBalance -A TandemRepeatAnnotator \
-A ClippingRankSumTest -A GCContent -A MappingQualityZero \
-A SpanningDeletions -A StrandOddsRatio -A AlleleBalanceBySample 
'''

def download_and_align(s3, bucket, fq1, fq2, sample, read_group_id, threads, ref, mem):
    fq1_local = "/ephemeral/" + os.path.basename(fq1)
    fq2_local = "/ephemeral/" + os.path.basename(fq2)
    
    logging.info("Downloading {} to {}".format(bucket + '/' + fq1, fq1_local))
    s3.Object(bucket, fq1).download_file(fq1_local)
    logging.info("Downloading {} to {}".format(bucket + '/' + fq2, fq2_local))
    s3.Object(bucket, fq2).download_file(fq2_local)
    
    cmd = align_cmd
    cmd = cmd.format(
        threads = threads,
        read_group = r"@RG\tID:{}\tSM:{}".format(read_group_id, sample),
        ref = ref,
        fq1 = fq1_local,
        fq2 = fq2_local,
        mem = mem,
        out = "/ephemeral/{read_group}_sorted.bam".format(read_group=read_group_id))

    logging.info("Running alignment: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    logging.info("Removing {} and {}".format(fq1_local, fq2_local))
    os.remove(fq1_local)
    os.remove(fq2_local)

    return "/ephemeral/{read_group}_sorted.bam".format(read_group=read_group_id)

def call_vars(bams, sample_name, ref, mem, gatk):
    hc_input = " -I " + " -I ".join(bams)
    hc_output = "/ephemeral/" + sample_name + ".g.vcf.gz"

    cmd = call_vars_cmd
    cmd = cmd.format(
        mem = mem,
        gatk = gatk,
        ref = ref,
        input = hc_input,
        out = hc_output)

    logging.info("Running variant calling: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    for bam in bams:
        logging.info("Removing {}".format(bam))
        os.remove(bam)
        logging.info("Removing {}".format(bam + ".bai"))
        os.remove(bam + ".bai")

    return hc_output

def index_bam(bam, threads):
    cmd = index_cmd
    cmd = cmd.format(threads=threads, bam=bam)

    logging.info("Indexing {}".format(bam))
    subprocess.check_call(cmd, shell=True)
    return

def process_args():
    parser = argparse.ArgumentParser(description="Process fastq files from a single sample to gvcf files. Initial files and final results are uploaded to S3")
    parser.add_argument("--threads", type=int, default=4, help="The number of alignment and indexing threads")
    parser.add_argument("--sort_mem", default="384M", help="Memory to use when sorting the alignment")
    parser.add_argument("--call_vars_mem", default="3g", help="Memory to use when calling variants")
    parser.add_argument("--gatk", default="/usr/local/bin/GenomeAnalysisTK.jar", help="The GATK .jar file")
    parser.add_argument("reference", help="The reference genome")
    parser.add_argument("access_key", help="AWS access key")
    parser.add_argument("secret_key", help="AWS secret key")
    parser.add_argument("upload_location", help="The S3 location for file upload")
    parser.add_argument("sample_name", help="The name of the sample")
    parser.add_argument("input_fastq", nargs='+', help="The input fastq files in s3")
    return parser.parse_args()

def main(args):
    log_format = "%(filename)s::%(funcName)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO,
                        format=log_format)

    if not args:
        args = process_args()

    session = boto3.session.Session(aws_access_key_id=args.access_key, aws_secret_access_key=args.secret_key)
    s3 = session.resource("s3")

    bams = []
    for fastq in args.input_fastq:
        bucket_end = fastq.find('/')
        bucket = fastq[:bucket_end]
        if bucket.startswith("s3://"):
            bucket = bucket[5:]
        fq1 = fastq[bucket_end + 1:]
        fq2 = fq1[:-15] + "2.filt.fastq.gz"
        read_group_id = os.path.basename(fq1)[:-16]

        # Download the fastq #
        next_bam = download_and_align(s3, bucket, fq1, fq2, args.sample_name, read_group_id, args.threads, args.reference, args.sort_mem)
        index_bam(next_bam, args.threads)
        bams.append(next_bam)

    # Call variants on the bam files #
    gvcf_local = call_vars(bams, args.sample_name, args.reference, args.call_vars_mem, args.gatk)
    gvcf_index = gvcf_local + ".tbi"

    # Upload the GVCF file #
    bucket_end = args.upload_location.find('/')
    bucket = args.upload_location[:bucket_end]
    if bucket.startswith("s3://"):
        bucket = bucket[5:]
    key = args.upload_location[bucket_end + 1:]
    logging.info("Uploading {} to {}".format(gvcf_local, bucket + '/' + key))
    s3.meta.client.upload_file(gvcf_local, bucket, key)
    logging.info("Uploading {} to {}".format(gvcf_index, bucket + '/' + key + ".tbi"))
    s3.meta.client.upload_file(gvcf_index, bucket, key + ".tbi")

    logging.info("Removing gvcf file: {}".format(gvcf_local))
    os.remove(gvcf_local)
    logging.info("Removing index file: {}".format(gvcf_index))
    os.remove(gvcf_index)

if __name__ == "__main__":
    main(None)
