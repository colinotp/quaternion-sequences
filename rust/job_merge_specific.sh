#! /bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=2000M
#SBATCH --account=def-cbright


# Script to run a SLURM job on a DRAC cluster
#
# This script goes merges a set of sorted files together

type=$1
filename=$2
output_dir=$3
n=$4

# sorting the files
start=`date +%s`
sort -m $output_dir/*.sorted > $filename.sorted
end=`date +%s`
if [[ $? -eq 0 ]]; then
    echo -e "Merging the file \"$(basename "$filename")\" took $((end - start)) seconds. \n\n" >> "./results/pairs/$type/find_$n/result.log"
fi
