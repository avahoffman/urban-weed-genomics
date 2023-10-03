#!/bin/bash

#SBATCH --time=15:00:00
#SBATCH --job-name=DS_STRUCT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=structure_ds.%j.out
#SBATCH --output=structure_ds.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4


# Written by: Ava Hoffman
# Date:       23 Feb 2023

export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/scratch_AH/structure_run

structure -m mainparams_ds1 -e extraparams_naive -o structure_out_DS1_naive
structure -m mainparams_ds2 -e extraparams_naive -o structure_out_DS2_naive
structure -m mainparams_ds3 -e extraparams_naive -o structure_out_DS3_naive
structure -m mainparams_ds4 -e extraparams_naive -o structure_out_DS4_naive
structure -m mainparams_ds5 -e extraparams_naive -o structure_out_DS5_naive
