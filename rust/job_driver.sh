#! /bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=4096M
#SBATCH --account=def-cbright

# Script to run a SLURM job on a DRAC cluster
#
# Simply calls the main driver script to compute QTS/PQS of a given length

./driver.sh qts 17 -dhs
