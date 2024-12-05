#!/bin/bash

#SBATCH --time=5:00:00
#SBATCH --job-name=PA_STRUCT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=structure_pa.%j.out
#SBATCH --output=structure_pa.%j.out
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4


# Written by: Ava Hoffman
# Date:       23 Feb 2023

export PATH=$PATH:/home/ahoffm31/scratch4-mavolio2/scratch_AH/structure_run

structure -m mainparams_pa1 -e extraparams_naive -o structure_out_PA1_naive_rep2
structure -m mainparams_pa2 -e extraparams_naive -o structure_out_PA2_naive_rep2
structure -m mainparams_pa3 -e extraparams_naive -o structure_out_PA3_naive_rep2
structure -m mainparams_pa4 -e extraparams_naive -o structure_out_PA4_naive_rep2
structure -m mainparams_pa5 -e extraparams_naive -o structure_out_PA5_naive_rep2

structure -m mainparams_pa1 -e extraparams_naive -o structure_out_PA1_naive_rep3
structure -m mainparams_pa2 -e extraparams_naive -o structure_out_PA2_naive_rep3
structure -m mainparams_pa3 -e extraparams_naive -o structure_out_PA3_naive_rep3
structure -m mainparams_pa4 -e extraparams_naive -o structure_out_PA4_naive_rep3
structure -m mainparams_pa5 -e extraparams_naive -o structure_out_PA5_naive_rep3

structure -m mainparams_pa1 -e extraparams_naive -o structure_out_PA1_naive_rep4
structure -m mainparams_pa2 -e extraparams_naive -o structure_out_PA2_naive_rep4
structure -m mainparams_pa3 -e extraparams_naive -o structure_out_PA3_naive_rep4
structure -m mainparams_pa4 -e extraparams_naive -o structure_out_PA4_naive_rep4
structure -m mainparams_pa5 -e extraparams_naive -o structure_out_PA5_naive_rep4

structure -m mainparams_pa1 -e extraparams_naive -o structure_out_PA1_naive_rep5
structure -m mainparams_pa2 -e extraparams_naive -o structure_out_PA2_naive_rep5
structure -m mainparams_pa3 -e extraparams_naive -o structure_out_PA3_naive_rep5
structure -m mainparams_pa4 -e extraparams_naive -o structure_out_PA4_naive_rep5
structure -m mainparams_pa5 -e extraparams_naive -o structure_out_PA5_naive_rep5