#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --job-name=popns
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
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
export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/code_AH/stacks/bin/
export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/code_AH/stacks-2.62/

while read p; do

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
  echo $p
  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
  popmap_file=$(ls /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/popmap*.txt)
  mkdir -p /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_10_pct/
  mkdir -p /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_20_pct/
  mkdir -p /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_25_pct/
  mkdir -p /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_30_pct/
  mkdir -p /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_35_pct/
  mkdir -p /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_40_pct/
  mkdir -p /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_45_pct/
  mkdir -p /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_50_pct/

  populations -P /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/ \
  -O /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_10_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 10 --min-maf 0.05
  # must be in a minimum of 10% of samples

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/ \
  -O /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_20_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 20 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/ \
  -O /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_25_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 25 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/ \
  -O /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_30_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 30 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/ \
  -O /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_35_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 35 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/ \
  -O /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_40_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 40 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/ \
  -O /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_45_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 45 --min-maf 0.05

  echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

  populations -P /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/ \
  -O /home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/"$p"/populations_50_pct/ \
  -M "$popmap_file" \
  -t 8 --structure --genepop --write-random-snp --verbose --vcf \
  --min-samples-overall 50 --min-maf 0.05

done </home/ahoffm31/scratch4-mavolio2/scratch_AH/catalog_tests/12-city_catalog_names.txt
