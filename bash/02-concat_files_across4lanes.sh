#!/bin/bash


# GIVE ABOUT 30 seconds FOR EACH SET OF 4 FILES
#SBATCH --time=0:15:00
#SBATCH --job-name=concat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=error_concat_acrosslanes.%j.out
#SBATCH --output=out_concat_acrosslanes.%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=1


# Written by: Ava Hoffman
# Date:       11 October 2021
# Purpose:    Concatenate raw sequencing data by sample that was run across 4 
#             lanes.

# Alternative approach can be found here:
# https://www.michaelchimenti.com/2015/08/concatenate-several-lanes-of-illumina-hiseq-reads-quickly/

# This approach takes the substring for each sub library and then concatenates 
# all 1-4 lanes sharing that substring. Assumes there are all 4 lanes (L001 
# through L004).
for i in *_L001_*.fastq.gz; do
b=${i%%_L001_*.fastq.gz}
# Read 1
cat "$b"_L001_R1_001.fastq.gz "$b"_L002_R1_001.fastq.gz "$b"_L003_R1_001.fastq.gz "$b"_L004_R1_001.fastq.gz > "$b"_R1.fastq.gz
# Read 2
cat "$b"_L001_R2_001.fastq.gz "$b"_L002_R2_001.fastq.gz "$b"_L003_R2_001.fastq.gz "$b"_L004_R2_001.fastq.gz > "$b"_R2.fastq.gz
done