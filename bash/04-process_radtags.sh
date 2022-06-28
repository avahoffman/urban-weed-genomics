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
# This creates a new directory for each sublibrary. This is where the output files and
# samples will be stored and is important so that the log file is not overwritten
# for each new sublibrary loop in the designated text file (i.e., 04-process_radtags_file_names.txt).
  out=/home-1/ahoffm31@jhu.edu/scratch/demux/sublibrary_"$p"
  mkdir -p $out
# Move into the new directory
  cd $out

#Process_radtags command
  process_radtags -1 ~/scratch/out_dir/"$p"_R1.1.fq.gz \
  -2 ~/scratch/out_dir/"$p"_R2.2.fq.gz \
  -i gzfastq -b ~/code/demux/"$p"_demux.txt \
  -o /home-1/ahoffm31@jhu.edu/scratch/demux -c -q \
  --inline_inline --renz_1 pstI --renz_2 mspI --rescue --disable_rad_check
done <~/code/04-process_radtags_file_names.txt
