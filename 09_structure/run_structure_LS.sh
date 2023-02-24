#!/bin/bash

#SBATCH --time=0:20:00
#SBATCH --job-name=LS_STRUCT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=structure_ls.%j.out
#SBATCH --output=structure_ls.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4


# Written by: Ava Hoffman
# Date:       23 Feb 2023

export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/scratch_AH/structure_run

structure -m mainparams_ls -e extraparams_naive -o structure_out_LS_naive
structure -m mainparams_ls -e extraparams_usepop -o structure_out_LS_popprior