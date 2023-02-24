#!/bin/bash

#SBATCH --time=0:50:00
#SBATCH --job-name=EC_STRUCT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=structure_ec.%j.out
#SBATCH --output=structure_ec.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4


# Written by: Ava Hoffman
# Date:       23 Feb 2023

export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/scratch_AH/structure_run

structure -m mainparams_ec -e extraparams_naive -o structure_out_EC_naive
structure -m mainparams_ec -e extraparams_usepop -o structure_out_EC_popprior