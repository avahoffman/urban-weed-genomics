#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --job-name=TO_STRUCT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=structure_to.%j.out
#SBATCH --output=structure_to.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4


# Written by: Ava Hoffman
# Date:       23 Feb 2023

export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/scratch_AH/structure_run

structure -m mainparams_to3_final -e extraparams_naive -o structure_out_TO3_naive_final