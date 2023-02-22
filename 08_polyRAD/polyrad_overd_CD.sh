#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --job-name=cd_OD
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=pr_cd.%j.out
#SBATCH --output=pr_cd.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4


# Written by: Ava Hoffman
# Date:       20 Feb 2023

module load r

Rscript polyRAD_overdispersion_CD.R