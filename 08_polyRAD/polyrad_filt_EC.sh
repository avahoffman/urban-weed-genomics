#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --job-name=ec_ft
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=pr_ec.%j.out
#SBATCH --output=pr_ec.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4


# Written by: Ava Hoffman
# Date:       20 Feb 2023

module load r

Rscript polyRAD_filter_EC.R