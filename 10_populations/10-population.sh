#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --job-name=popns
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=populations.%j.out
#SBATCH --output=populations.%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=2


# Written by: Ava Hoffman
# Date:       24 October 2021
# Purpose:    Running populations
# populations -P dir [-O dir] [-M popmap] (filters) [--fstats] [-k [--sigma=150000] [--bootstrap [-N 100]]] (output formats)

module load gcc
export PATH=$PATH:/home-1/ahoffm31@jhu.edu/code/bin

populations -P ~/scratch/stacks/stacks_21_11_09 -M ~/code/08-tsv2bam_popmap.txt -t 18 \
--structure --genepop --write-single-snp --verbose --vcf -r 0.1
