#!/bin/bash

#SBATCH --time=5:00:00
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

structure -m mainparams_ls1 -e extraparams_naive -o structure_out_LS1_naive
#tructure -m mainparams_ls -e extraparams_usepop -o structure_out_LS_popprior
structure -m mainparams_ls2 -e extraparams_naive -o structure_out_LS2_naive
structure -m mainparams_ls3 -e extraparams_naive -o structure_out_LS3_naive
structure -m mainparams_ls4 -e extraparams_naive -o structure_out_LS4_naive
structure -m mainparams_ls5 -e extraparams_naive -o structure_out_LS5_naive