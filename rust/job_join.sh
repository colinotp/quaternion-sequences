#! /bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=4000M
#SBATCH --account=def-cbright

# Script to run a SLURM job on a DRAC cluster
#
# This script runs the part of the algorithm that goes
# through the sorted auto and cross correlation values
# to find valid QTS, and then computes the corresponding PQS

type=$1
n=$2

./join_pairs.sh $type $n

