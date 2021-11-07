#!/bin/bash

# Need 10 mins per sample
#SBATCH --time=05:00:00
#SBATCH --job-name=uStcks
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=ustacks1.%j.out
#SBATCH --output=ustacks1.%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=1


# Written by: Ava Hoffman
# Date:       21 October 2021
# Purpose:    Running ustacks on each of the samples specified, building loci de novo in each sample
# ustacks -f file_path  -i id -o ~/scratch/stacks/ustacks [-M max_dist] [-m min_cov] [-p num_threads]

module load gcc
export PATH=$PATH:/home-1/ahoffm31@jhu.edu/code/bin

while read -r f1 f2; do
  ustacks -f ~/scratch/demux/DS/"$f1".1.fq.gz -i "$f2" \
  -o ~/scratch/stacks/ustacks -M 2 -m 3 -N 4 -p 16 \
  --max_locus_stacks 3 --deleverage
done < <(paste ~/code/05-ustacks_samples_1.txt ~/code/05-ustacks_id_1.txt)