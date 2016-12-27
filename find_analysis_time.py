#!/usr/bin/env python3

import argparse
import logging
import re
import datetime

def process_args():
    parser = argparse.ArgumentParser(description="Find the compute time required to process the 1000 genomes dataset")
    parser.add_argument("sequence_index")
    parser.add_argument("log_files", nargs='+')
    return parser.parse_args()

def datetime_from_line(line):
    line = line.split()[:2]
    return datetime.datetime(*[int(x) for x in line[0].split('-')] + [int(x) for x in re.split(r"[:,]", line[1])])

def main(args):
    if not args:
        args = process_args()

    samples = {}
    bam = {}
    with open(args.sequence_index) as f:
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
            base_count = line[header.index("BASE_COUNT")]


            if not sample in samples:
                samples[sample] = []
            samples[sample].append(run_id)
            bam[run_id] = int(base_count)

    alignment_runtimes = []
    dedup_runtimes = []
    vc_runtimes = {}
    cat_runtimes = {}

    for log_file in args.log_files:
        with open(log_file, 'r') as f:
            run_datetime = None
            run_id = None
            run_id2 = None
            sample_id = None
            cat_sample_id = None
            for line in f:
                if "Running alignment:" in line:
                    run_datetime = datetime_from_line(line)
                    cmd = f.readline()
                    run_id = re.search(r'\\tID:([SE]RR[0-9]+)\\tSM', cmd)
                    if run_id:
                        run_id = run_id.group(1)

                elif "Running duplicate removal:" in line:
                    run_datetime = datetime_from_line(line)
                    cmd = f.readline()
                    run_id2 = re.search(r'/ephemeral/([SE]RR[0-9]+)_dup.bam', cmd)
                    if run_id2:
                        run_id2 = run_id2.group(1)

                elif "Running variant calling:" in line:
                    run_datetime = datetime_from_line(line)
                    cmd = f.readline()
                    sample_id = re.search(" /ephemeral/(HG[0-9]+)", cmd)
                    if sample_id:
                        sample_id = sample_id.group(1)

                elif "Running cat variants:" in line:
                    run_datetime = datetime_from_line(line)
                    cmd = f.readline()
                    cat_sample_id = re.search(" -out /ephemeral/(HG[0-9]+)", cmd)
                    if cat_sample_id:
                        cat_sample_id = cat_sample_id.group(1)

                elif "[INFO]" in line and run_datetime and run_id:
                    cur_datetime = datetime_from_line(line)
                    alignment_runtimes.append((bam[run_id], (cur_datetime - run_datetime).total_seconds()))
                    run_id = None
                    run_datetime = None

                elif "[INFO]" in line and run_datetime and run_id2:
                    cur_datetime = datetime_from_line(line)
                    dedup_runtimes.append((bam[run_id2], (cur_datetime - run_datetime).total_seconds()))
                    run_id2 = None
                    run_datetime = None

                elif "[INFO]" in line and run_datetime and sample_id:
                    cur_datetime = datetime_from_line(line)
                    if not sample_id in vc_runtimes:
                        vc_runtimes[sample_id] = []
                    vc_runtimes[sample_id].append((cur_datetime - run_datetime).total_seconds())
                    sample_id = None
                    run_datetime = None

                elif "[INFO]" in line and run_datetime and cat_sample_id:
                    cur_datetime = datetime_from_line(line)
                    if not cat_sample_id in cat_runtimes:
                        cat_runtimes[cat_sample_id] = []
                    cat_runtimes[cat_sample_id] = (cur_datetime - run_datetime).total_seconds()
                    cat_sample_id = None
                    run_datetime = None
        
                elif "Traceback (most recent call last):" in line:
                    run_datetime = None
                    run_id = None
                    sample_id = None
                    cat_sample_id = None

    try:
        aligned_bases, alignment_time = list(map(sum, zip(*alignment_runtimes)))
    except TypeError:
        print(str(alignment_runtimes))
        raise
    rate = aligned_bases / alignment_time
    all_bases = sum(bam.values())
    print("Aligned {} bases in {} seconds or {} hours".format(aligned_bases, alignment_time, alignment_time / 3600))
    print("All {} bases would take {} hours".format(all_bases, all_bases / rate / 3600))
   
    if dedup_runtimes:
        try:
            dedup_bases, dedup_time = list(map(sum, zip(*dedup_runtimes)))
        except TypeError:
            print(str(dedup_runtimes))
            raise
        rate = dedup_bases / dedup_time
        all_bases = sum(bam.values())
        print("Deduped {} bases in {} seconds or {} hours".format(dedup_bases, dedup_time, dedup_time / 3600))
        print("All {} bases would take {} hours".format(all_bases, all_bases / rate / 3600))

    total_vc_runtimes = []
    for sample, runtimes in vc_runtimes.items():
        if len(runtimes) == 25 or len(runtimes) == 1: # All chromosomes present or a single Sentieon run
            total_vc_runtimes.append((sum([bam[x] for x in samples[sample]]), sum(runtimes)))
    bases, time = list(map(sum, zip(*total_vc_runtimes)))
    rate = bases / time
    print("Called {} bases in {} seconds or {} hours".format(bases, time, time / 3600))
    print("All {} bases would take {} hours".format(all_bases, all_bases / rate / 3600))

    total_cat_runtimes = []
    for sample, runtime in cat_runtimes.items():
        total_cat_runtimes.append((sum([bam[x] for x in samples[sample]]), runtime))
    if total_cat_runtimes:
        bases, time = list(map(sum, zip(*total_cat_runtimes)))
        rate = bases / time
        print("Concatenated {} bases in {} seconds or {} hours".format(bases, time, time / 3600))
        print("All {} bases would take {} hours".format(all_bases, all_bases / rate / 3600))

if __name__ == "__main__":
    main(None)
