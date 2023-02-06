#!/bin/bash

# Estimated 4 minutes (max) per sample @ tasks =3 and -p 24 BUT this increases
# With number of samples as the catalog gets more complex.
# 45 minutes for 30 samples
#SBATCH --time=4:00:00
#SBATCH --job-name=cStacks
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jenny.cocciardi@jhu.edu
#SBATCH --error=cstacks.%j.out
#SBATCH --output=cstacks.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=6


# Written by: Jenny Cocciardi
# Date:       13 July 2022
# Purpose:    Running cstacks to create a set of consensus loci, merging alleles together
# cstacks -P in_path -M popmap -o output [-n num_mismatches] [-p num_threads]

module load gcc
export PATH=$PATH:/home/jcoccia1/data_mavolio2/code/test/stacks-2.60/local/bin

cstacks -P ~/scr4_mavolio2/scratch_JC/stacks/CD/catalog_7_13_22 \
  -M ~/scr4_mavolio2/scratch_JC/stacks/CD/catalog_7_13_22/cstacks_popmap_CD.txt \
  -n 8 -p 24
