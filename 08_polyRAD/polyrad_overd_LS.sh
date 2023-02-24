#!/bin/bash

#SBATCH --time=0:05:00
#SBATCH --job-name=ls_OD
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=pr_ls.%j.out
#SBATCH --output=pr_ls.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2


# Written by: Ava Hoffman
# Date:       20 Feb 2023

module load r

Rscript polyRAD_overdispersion_LS.R