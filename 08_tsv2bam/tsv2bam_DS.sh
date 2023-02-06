#!/bin/bash

# 19 mins for 107 samples
# 32 mins for 221 samples
#SBATCH --time=00:45:00
#SBATCH --job-name=t2bam
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jenny.cocciardi@jhu.edu
#SBATCH --error=t2bam.%j.out
#SBATCH --output=t2bam.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=5


# Written by: Ava Hoffman
# Date:       24 October 2021
# Purpose:    Running tsv2bam to transpose data so that it is oriented by locus, instead of by sample
# tsv2bam -P ./stacks/ -M ./popmap -R ./samples -t 8

module load gcc
export PATH=$PATH:/home/jcoccia1/data_mavolio2/code/test/stacks-2.60/local/bin

tsv2bam -P ~/scr4_mavolio2/scratch_JC/stacks/DS/catalog_7_7_22 \
 -M ~/scr4_mavolio2/scratch_JC/stacks/DS/catalog_7_7_22/popmap_DS.txt \
 -R ~/scr4_mavolio2/scratch_JC/demux/DS -t 18
