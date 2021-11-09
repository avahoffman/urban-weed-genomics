#!/bin/bash

# ~2-2.25 minutes for each sample, with -p = 8 and ntasks = 3
# But will depend on the size of the catalog.
# Can be done piecemeal.
#SBATCH --time=09:00:00
#SBATCH --job-name=sStcks
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=sstacks.%j.out
#SBATCH --output=sstacks.%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=2
# Good to have a lot of memory on this one
# when looping, do 1 hr, with 10 tasks per node


# Written by: Ava Hoffman
# Date:       22 October 2021
# Purpose:    All samples in the populations are matched against the catalog
# sstacks -c catalog_path -s sample_path [-s sample_path ...] -o path [-p n_threads]

module load gcc
export PATH=$PATH:/home-1/ahoffm31@jhu.edu/code/bin

# N=4
# (
# while read p; do
#   ((i=i%N)); ((i++==0)) && wait
#   sstacks -c ~/scratch/stacks/ustacks -s ~/scratch/stacks/ustacks/"$p" \
#   -o ~/scratch/stacks/sstacks -p 8 &
# done <~/code/sstacks_samples.txt
# )

while read p; do
sstacks -c ~/scratch/stacks/sstacks_21_11_09 -s ~/scratch/stacks/sstacks_21_11_09/"$p" \
  -o ~/scratch/stacks/sstacks_21_11_09 -p 8
done <~/code/07-sstacks_samples.txt