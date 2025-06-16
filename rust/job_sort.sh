#! /bin/bash
#SBATCH --time=30:00:00
#SBATCH --mem=4000M
#SBATCH --account=def-cbright

# Script to run a SLURM job on a DRAC cluster
#
# This script runs sorts the .pair files generated in the first part of the algorithm.
# Example usage: sbatch job_sort.sh n
# where n is the desired length

n=$1

# sorting the files
start2=`date +%s`
./sortpairs.sh qts $n
end2=`date +%s`
echo Sorting the files took `expr $end2 - $start2` seconds.
