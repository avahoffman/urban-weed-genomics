#!/bin/bash

# One hour per sublibrary (?)
#SBATCH --time=14:00:00
#SBATCH --job-name=cleanmx
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=process_radtags.%j.out
#SBATCH --output=process_radtags.%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=1


# Written by: Ava Hoffman
# Date:       11 October 2021
# Purpose:    Run process_radtags to demultiplex and clean for quality

module load gcc
export PATH=$PATH:/home-1/ahoffm31@jhu.edu/code/bin

while read p; do
  process_radtags -1 ~/scratch/out_dir/"$p"_R1.1.fq.gz \
  -2 ~/scratch/out_dir/"$p"_R2.2.fq.gz \
  -i gzfastq -b ~/code/demux/"$p"_demux.txt \
  -o /home-1/ahoffm31@jhu.edu/scratch/demux -c -q \
  --inline_inline --renz_1 pstI --renz_2 mspI --rescue --disable_rad_check
done <~/code/process_radtags_file_names.txt