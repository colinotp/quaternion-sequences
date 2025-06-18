#! /bin/bash
#SBATCH --time=30:00:00
#SBATCH --mem=4000M
#SBATCH --account=def-cbright

# Script to run a SLURM job on a DRAC cluster
#
# This script sorts the .pair files generated in the first part of the algorithm,
# using a combination of merge sort and bash quick sort for handling of large files.

# Example usage: sbatch job_sort_split.sh n
# where n is the desired length

type=$1
n=$2
./sortpairs_split.sh $type $n
