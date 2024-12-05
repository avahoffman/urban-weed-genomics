#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH --job-name=CD_STRUCT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=structure_cd.%j.out
#SBATCH --output=structure_cd.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4


# Written by: Ava Hoffman
# Date:       23 Feb 2023

export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/scratch_AH/structure_run

structure -m mainparams_cd3_final -e extraparams_naive -o structure_out_CD3_naive_final