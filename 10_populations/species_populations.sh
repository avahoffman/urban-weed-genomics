#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --job-name=popns
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jenny.cocciardi@jhu.edu
#SBATCH --error=populations.%j.out
#SBATCH --output=populations.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=2


# Written by: Ava Hoffman
# Date:       27 September 2022
# Purpose:    Running populations across various 
# populations -P dir [-O dir] [-M popmap] (filters) [--fstats] [-k [--sigma=150000] [--bootstrap [-N 100]]] (output formats)

module load gcc
export PATH=$PATH:/home/jcoccia1/data_mavolio2/code/test/stacks-2.60/local/bin

while read p; do

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
  echo $p
  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
  popmap_file=$(ls /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/popmap*.txt)
  mkdir -p /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_10_pct/
  mkdir -p /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_20_pct/
  mkdir -p /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_25_pct/
  mkdir -p /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_30_pct/
  mkdir -p /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_35_pct/
  mkdir -p /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_40_pct/
  mkdir -p /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_45_pct/
  mkdir -p /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_50_pct/

  populations -P /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/ \
  -O /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_10_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 10 --min-maf 0.05
  # must be in a minimum of 10% of samples

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/ \
  -O /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_20_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 20 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/ \
  -O /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_25_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 25 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/ \
  -O /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_30_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 30 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/ \
  -O /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_35_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 35 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/ \
  -O /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_40_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 40 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/ \
  -O /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_45_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 45 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/ \
  -O /home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/"$p"/populations_50_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 50 --min-maf 0.05

done </home/jcoccia1/scr4_mavolio2/scratch_JC/metapop_catalogs/species_catalog_names.txt
