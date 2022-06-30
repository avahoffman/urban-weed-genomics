#!/bin/bash

# Need  ~15 mins per sample
#SBATCH --time=3:00:00
#SBATCH --job-name=denovo_map
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jcoccia1@jhu.edu
#SBATCH --error=denovo_map_CD.%j.out
#SBATCH --output=denovo_map_CD.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=2


# Written by: Jenny Cocciardi
# Date:       21 March 2022
# Purpose:    Running denovo_map.pl to execute the Stacks pipeline on a subset of samples for \
#             parameter selection - doing this for CD species
# denovo_map.pl --samples dir --popmap path -o dir [--paired [--rm-pcr-duplicates]] \
#       (assembly options) (filtering options) [-X prog:"opts" ...]

module load gcc
module load bioperl
export PATH=$PATH:/home/jcoccia1/data_mavolio2/code/test/stacks-2.60/local/bin

denovo_map.pl -M 2 -o ~/scr4_mavolio2/scratch_JC/04b-demux/CD/05a-tests_denovo/stacks.M2 \
  --popmap ~/scr4_mavolio2/scratch_JC/04b-demux/CD/05a-tests_denovo/stacks.M2/r80/r80_2/popmap_test_samples.txt \
  --samples ~/scr4_mavolio2/scratch_JC/04b-demux/CD --paired --min-samples-per-pop 0.80

