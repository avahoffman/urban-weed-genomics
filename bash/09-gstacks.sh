#!/bin/bash

# 57 mins for 107 samples @ 2 tasks per node
#SBATCH --time=02:20:00
#SBATCH --job-name=gStcks
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=gstacks.%j.out
#SBATCH --output=gstacks.%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=2


# Written by: Ava Hoffman
# Date:       24 October 2021
# Purpose:    Running gstacks
# gstacks -P stacks_dir -M popmap

module load gcc
export PATH=$PATH:/home-1/ahoffm31@jhu.edu/code/bin

gstacks -P ~/scratch/stacks/stacks_21_11_09 -M ~/code/08-tsv2bam_popmap.txt -t 18