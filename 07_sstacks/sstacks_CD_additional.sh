#!/bin/bash

# ~2-2.25 minutes for each sample, with -p = 8 and ntasks = 3
 # But will depend on the size of the catalog.
# Can be done piecemeal.
#SBATCH --time=2:30:00
#SBATCH --job-name=sstacks
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jenny.cocciardi@jhu.edu
#SBATCH --error=sstacks.%j.out
#SBATCH --output=sstacks.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=4
# Good to have a lot of memory on this one
# when looping, do 1 hr, with 10 tasks per node


# Written by: Ava Hoffman
# Date:       22 October 2021
# Purpose:    All samples in the populations are matched against the catalog
# sstacks -c catalog_path -s sample_path [-s sample_path ...] -o path [-p n_threads]

module load gcc
export PATH=$PATH:/home/jcoccia1/data_mavolio2/code/test/stacks-2.60/local/bin

while read p; do
sstacks -c ~/scr4_mavolio2/scratch_JC/stacks/CD/CD_metapop_catalog \
  -s ~/scr4_mavolio2/scratch_JC/stacks/CD/CD_metapop_catalog/"$p" \
  -o ~/scr4_mavolio2/scratch_JC/stacks/CD/CD_metapop_catalog \
  -p 10
done <~/scr4_mavolio2/scratch_JC/stacks/CD/CD_metapop_catalog/sstacks_samples_CD_additional.txt
