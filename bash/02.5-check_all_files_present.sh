#!/bin/bash


# Need at max 1 minute
#SBATCH --time=00:01:00
#SBATCH --job-name=filecheck
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=error_filecheck.%j.out
#SBATCH --output=out_filecheck.%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=1


# Written by: Ava Hoffman
# Date:       9 February 2022
# Purpose:    Check that all files made it via Aspera transfer and were 
#             concatenated as expected.
# Useage:     Run from the directory where concatenated files are stored (e.g.,
#             the scratch directory). Do not edit 02.5-all_possible_clonefilter_filenames.txt 
#             unless checking that a subset of files are present.

while read p; do
  if [ -f "$p"_R1.fastq.gz ];
    then
      continue
    else
	    echo "$p"_R1.fastq.gz file does not exist.
  fi
  if [ -f "$p"_R2.fastq.gz ];
    then
      continue
    else
	    echo "$p"_R2.fastq.gz file does not exist.
  fi
done <~/code/02.5-all_possible_clonefilter_filenames.txt