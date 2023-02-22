#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --job-name=ls_polyRAD
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=pr_ls.%j.out
#SBATCH --output=pr_ls.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


# Written by: Ava Hoffman
# Date:       20 Feb 2023

module load r

Rscript make_polyRAD_LS.R