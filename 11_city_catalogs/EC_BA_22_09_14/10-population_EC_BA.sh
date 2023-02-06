#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --job-name=EC_BA_popns
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=populations.%j.out
#SBATCH --output=populations.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=2


# Written by: Ava Hoffman
# Date:       14 September 2022
# Purpose:    Running populations
# populations -P dir [-O dir] [-M popmap] (filters) [--fstats] [-k [--sigma=150000] [--bootstrap [-N 100]]] (output formats)

module load gcc
export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/code_AH/stacks/bin/
export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/code_AH/stacks-2.62/

populations -P /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/EC_BA_22_09_14/ \
  -M /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/EC_BA_22_09_14/popmap_EC_BA.txt \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf
