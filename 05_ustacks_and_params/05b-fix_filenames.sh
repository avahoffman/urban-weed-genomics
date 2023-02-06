#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --job-name=filefx
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=filefix.%j.out
#SBATCH --output=filefix.%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
# 24 cores = a single node
#SBATCH --ntasks-per-node=2


# Written by: Ava Hoffman
# Date:       24 October 2021
# Purpose:    Fix filenames (stacks adds an extra "1" at some point. The following:
# DS.MN.L02-DS.M.3.1.alleles.tsv.gz
# DS.MN.L03-DS.U.2.1.tags.tsv.gz
# DS.MN.L09-DS.U.1.1.snps.tsv.gz
# become:
# DS.MN.L02-DS.M.3.alleles.tsv.gz
# DS.MN.L03-DS.U.2.tags.tsv.gz
# DS.MN.L09-DS.U.1.snps.tsv.gz

module load gcc
export PATH=$PATH:/home-1/ahoffm31@jhu.edu/code/bin

########3
for i in *; do mv -- "$i" "${i//1.1./1.}"; done
for i in *; do mv -- "$i" "${i//2.1./2.}"; done
for i in *; do mv -- "$i" "${i//3.1./3.}"; done
for i in *; do mv -- "$i" "${i//4.1./4.}"; done
for i in *; do mv -- "$i" "${i//5.1./5.}"; done
for i in *; do mv -- "$i" "${i//6.1./6.}"; done
for i in *; do mv -- "$i" "${i//7.1./7.}"; done