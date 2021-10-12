#!/bin/bash

# ABOUT 1.5 hour for each set of reads per library
#SBATCH --time=6:00:00
#SBATCH --job-name=clonefilter
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=clone_filter.%j.out
#SBATCH --output=clone_filter.%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=24


# Written by: Ava Hoffman
# Date:       11 October 2021
# Purpose:    Run clone_filter to remove pcr duplicates

module load gcc
export PATH=$PATH:/home-1/ahoffm31@jhu.edu/code/bin

# clone_filter -1 AMH_macro_1_2_12px_S2_R1.fastq.gz -2 AMH_macro_1_2_12px_S2_R2.fastq.gz -i gzfastq -o /home-1/ahoffm31@jhu.edu/scratch/out_dir --inline_inline --oligo_len_1 4 --oligo_len_2 4
# clone_filter -1 AMH_macro_8_1_8px_S99_R1.fastq.gz -2 AMH_macro_8_1_8px_S99_R2.fastq.gz -i gzfastq -o /home-1/ahoffm31@jhu.edu/scratch/out_dir --inline_inline --oligo_len_1 4 --oligo_len_2 4
# clone_filter -1 AMH_macro_8_2_8px_S100_R1.fastq.gz -2 AMH_macro_8_2_8px_S100_R2.fastq.gz -i gzfastq -o /home-1/ahoffm31@jhu.edu/scratch/out_dir --inline_inline --oligo_len_1 4 --oligo_len_2 4
# clone_filter -1 AMH_macro_8_3_8px_S101_R1.fastq.gz -2 AMH_macro_8_3_8px_S101_S2_R2.fastq.gz -i gzfastq -o /home-1/ahoffm31@jhu.edu/scratch/out_dir --inline_inline --oligo_len_1 4 --oligo_len_2 4

while read p; do
  clone_filter -1 "$p"_R1.fastq.gz -2 "$p"_R2.fastq.gz \
  -i gzfastq -o /home-1/ahoffm31@jhu.edu/scratch/out_dir \
  --inline_inline --oligo_len_1 4 --oligo_len_2 4
done <~/code/clone_filter_file_names.txt