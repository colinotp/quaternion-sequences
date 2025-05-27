#! /bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=2000M
#SBATCH --account=def-cbright

# Script to run a SLURM job on a DRAC cluster
#
# Simply calls the main driver script to compute QTS/PQS of a given length

./driver.sh 17
