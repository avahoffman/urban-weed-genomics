#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --job-name=ds_OD
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=pr_ds.%j.out
#SBATCH --output=pr_ds.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4


# Written by: Ava Hoffman
# Date:       20 Feb 2023

module load r

Rscript polyRAD_overdispersion_DS.R