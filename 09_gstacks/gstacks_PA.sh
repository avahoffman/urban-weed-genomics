#!/bin/bash

# 57 mins for 107 samples @ 2 tasks per node
#SBATCH --time=04:00:00
#SBATCH --job-name=gstacks
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jenny.cocciardi@jhu.edu
#SBATCH --error=gstacks.%j.out
#SBATCH --output=gstacks.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=3


# Written by: Jenny Cocciardi
# Date:       July 7 2022
# Purpose:    Running gstacks
# gstacks -P stacks_dir -M popmap

module load gcc
export PATH=$PATH:/home/jcoccia1/data_mavolio2/code/test/stacks-2.60/local/bin

gstacks -P ~/scr4_mavolio2/scratch_JC/stacks/PA/catalog_7_13_22 \
 -M ~/scr4_mavolio2/scratch_JC/stacks/PA/catalog_7_13_22/popmap_PA.txt -t 18
