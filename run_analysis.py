#!/usr/bin/env python3

import argparse
import subprocess
import boto3
import logging
import os.path
import time
import pickle
from retrying import retry

analysis_cmd = 'qsub -e {log} -o {out} -l h="node*",ephemeral={size},total_mem={mem} general.sh {fastq_to_gvcf} --bam_key {bam_key} {ref} {access_key} {secret_key} {dest} {sample} {input_fastq}'

@retry(wait_fixed=2000, stop_max_attempt_number=5)
def check_n_waiting_jobs(max_waiting_jobs):
    while True:
        out = subprocess.check_output("qstat", shell=True, universal_newlines=True)
        if not out: # No running jobs
            return
        out = out.split('\n')
        waiting = 0
        for job in out:
            if ' qw ' in job:
                waiting += 1
        if waiting < max_waiting_jobs:
            return

def process_args():
    parser = argparse.ArgumentParser(description="Schedule processing of 1000 genomes fastq to gVCF files")
    parser.add_argument("--fastq_index", default="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.sequence.index", help="The index of fastq files used in phase3")
    parser.add_argument("--onekg_bucket", default="1000genomes", help="The 1000 genomes bucket")
    parser.add_argument("--n_to_run", type=int, help="The number of samples to run [all]")
    parser.add_argument("--log_dir", default="/data/Logs/", help="A directory to store log files")
    parser.add_argument("--reference", default="/data/Reference/hs37d5.fa", help="The human reference genome")
    parser.add_argument("--sleep", default=5.0, type=float, help="The amount of time to sleep inbetween queueing jobs")
    parser.add_argument("--max_waiting_jobs", default=10, type=int, help="The maximum number of waiting jobs in the queue")
    parser.add_argument("--destination_key", default="1000genomes/gVCF/{sample}/{sample}.g.vcf.gz", help="The S3 destination key")
    parser.add_argument("--bam_key", default="1000genomes/BAM/{sample}/{run}.bam", help="The S3 destination for temporary BAM files")
    parser.add_argument("--s3_keys_cache", default="/home/onekg/s3_paths.p", help="A pickle of the fastq keys in s3")
    parser.add_argument("destination_bucket", help="The destination bucket")
    parser.add_argument("access_key", help="AWS access key")
    parser.add_argument("secret_key", help="AWS secret key")
    return parser.parse_args()

def main(args):
    log_format = "%(filename)s::%(funcName)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO,
                        format=log_format)

    if not args:
        args = process_args()

    aws_session = boto3.session.Session(aws_access_key_id=args.access_key, aws_secret_access_key=args.secret_key)
    s3 = aws_session.resource("s3")

    # Find all of the phase3 fastq in the 1000 genomes bucket #
    if os.path.isfile(args.s3_keys_cache):
        s3_paths = pickle.load(open(args.s3_keys_cache, "rb"))
    else:
        s3_paths = {}
        sample_bucket = s3.Bucket(args.onekg_bucket)
        for s3_obj in sample_bucket.objects.filter(Prefix="phase3/data/"):
            if s3_obj.key.endswith("_1.filt.fastq.gz"):
                fq_name = s3_obj.key.split('/')[-1][:-16]
                s3_paths[fq_name] = (s3_obj.key, s3_obj.size)
        pickle.dump(s3_paths, open(args.s3_keys_cache, "wb"))

    # Download the sequence index #
    if not os.path.isfile(os.path.basename(args.fastq_index)):
        cmd = "wget {}".format(args.fastq_index)
        subprocess.check_call(cmd, shell=True)
    
    # Get all fastq for a single sample #
    samples = {}
    with open(os.path.basename(args.fastq_index)) as f:
        header = f.readline().rstrip().split('\t')
        for line in f:
            line = line.rstrip().split('\t')
            if line[header.index("ANALYSIS_GROUP")] != "low coverage":
                continue
            if line[header.index("INSTRUMENT_PLATFORM")] != "ILLUMINA":
                continue
            if line[header.index("WITHDRAWN")] != "0":
                continue
            if line[header.index("LIBRARY_LAYOUT")] != "PAIRED":
                continue
            if not line[header.index("FASTQ_FILE")].endswith("_1.filt.fastq.gz"):
                continue
            
            sample = line[header.index("SAMPLE_NAME")]
            run_id = line[header.index("RUN_ID")]
            
            if not sample in samples:
                samples[sample] = []
            samples[sample].append(s3_paths[run_id])

    # Find all of the finsihed samples #
    finished_samples = set()
    dest_bucket = s3.Bucket(args.destination_bucket)
    for s3_obj in dest_bucket.objects.filter(Prefix=args.destination_key[:args.destination_key.index('{')]):
        if s3_obj.key.endswith(".g.vcf.gz"):
            finished_samples.add(s3_obj.key)            

    # Submit the individual jobs to the scheduler #
    n_run = 0
    for sample in sorted(samples.keys()):
        sample_key = args.destination_key.format(sample=sample)
        if sample_key in finished_samples:
            logging.info("Sample key {} is already present in the destination bucket {}. Skipping...".format(sample_key, args.destination_bucket))
            continue

        check_n_waiting_jobs(args.max_waiting_jobs)

        fastq, sizes = tuple(zip(*samples[sample]))
        
        cmd = analysis_cmd.format(
            log = args.log_dir + "/analysis_{}.log".format(n_run),
            out = args.log_dir + "/analysis_{}.out".format(n_run),
            size = sum(sizes) * 2,
            mem = 8 * 1024 * 1024 * 1024, # Alignment ~6 GB, sorting ~1 GB, duplicates ?
            fastq_to_gvcf = '/data/sample_fastq_to_gvcf.py',
            bam_key = args.bam_key,
            ref = args.reference,
            access_key = args.access_key,
            secret_key = args.secret_key,
            dest = args.destination_bucket + '/' + sample_key,
            sample = sample,
            input_fastq = ' '.join([args.onekg_bucket + '/' + x for x in fastq]))

        logging.info("Running:\n{}".format(cmd))
        subprocess.check_call(cmd, shell=True)
        time.sleep(args.sleep)

        n_run += 1
        if hasattr(args, "n_to_run"):
            if n_run >= args.n_to_run:
                break

if __name__ == "__main__":
    main(None)
