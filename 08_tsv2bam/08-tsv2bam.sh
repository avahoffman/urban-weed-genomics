#!/bin/bash

# 19 mins for 107 samples
# 32 mins for 221 samples
#SBATCH --time=00:45:00
#SBATCH --job-name=t2bam
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=t2bam.%j.out
#SBATCH --output=t2bam.%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=2


# Written by: Ava Hoffman
# Date:       24 October 2021
# Purpose:    Running tsv2bam to transpose data so that it is oriented by locus, instead of by sample
# tsv2bam -P ./stacks/ -M ./popmap -R ./samples -t 8

module load gcc
export PATH=$PATH:/home-1/ahoffm31@jhu.edu/code/bin

tsv2bam -P ~/scratch/stacks/stacks_21_11_09 -M ~/code/08-tsv2bam_popmap.txt \
-R ~/scratch/demux/DS -t 18