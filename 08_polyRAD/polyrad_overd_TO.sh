#!/bin/bash

#SBATCH --time=0:20:00
#SBATCH --job-name=to_OD
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=pr_to.%j.out
#SBATCH --output=pr_to.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4


# Written by: Ava Hoffman
# Date:       20 Feb 2023

module load r

Rscript polyRAD_overdispersion_TO.R