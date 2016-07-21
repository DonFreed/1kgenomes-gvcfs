# 1kgenomes-gvcfs

### Overview

This repository contains scripts for analyzing low-coverage fastq files using in Phase 3 of the 1000 Genomes project to produce gVCF files. The scripts are designed to run on Amazon EC2 using [StarCluster](http://star.mit.edu/cluster/). The scripts are run on a custom AMI which is described below. Raw fastq files from each sample are aligned with bwa, duplicates are marked with samblaster, alignments are sorted and indexed using sambamba and variants are called using the GATK. The `/data` directory contains an NFS shared EBS volume with the indexed reference genome ([hs37d5.fa](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/README_human_reference_20110707)) and space for log files.

### AMI Details

The AMI contains the following software:
* Python3
* * retrying
* * awscli
* * boto3
* mdadm (software RAID)
* [bwa](https://github.com/lh3/bwa) version 0.7.15
* [samblaster](https://github.com/GregoryFaust/samblaster) version 0.1.22
* [sambamba](http://lomereiter.github.io/sambamba/) version 0.6.3
* The [GATK](https://software.broadinstitute.org/gatk/) jar file at /home/ec2-user/ version 3.5.0
