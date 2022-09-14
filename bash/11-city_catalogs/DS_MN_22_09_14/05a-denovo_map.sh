#!/bin/bash

# Need  ~15 mins per sample
#SBATCH --time=12:00:00
#SBATCH --job-name=dm_DS_MN
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=denovo_map_DS_MN.%j.out
#SBATCH --output=denovo_map_DS_MN.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=8


# Written by: Ava Hoffman
# Date:       14 September 2022
# Purpose:    Running denovo_map.pl to execute the Stacks pipeline on a subset of samples for \
#             parameter selection - doing this for DS species
# denovo_map.pl --samples dir --popmap path -o dir [--paired [--rm-pcr-duplicates]] \
#       (assembly options) (filtering options) [-X prog:"opts" ...]

module load gcc
module load bioperl
export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/code_AH/stacks/bin/
export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/code_AH/stacks-2.62/

echo "Creating a by-city catalog for Digitaria (DS) in Minneapolis."

denovo_map.pl -M 10 -T 8 -o /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/DS_MN_22_09_14 \
  --popmap /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/DS_MN_22_09_14/popmap_DS_MN.txt \
  --samples /home/ahoffm31/scratch4-mavolio2/scratch_AH/demux/DS --paired --time-components

