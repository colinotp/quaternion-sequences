#! /bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=2000M
#SBATCH --account=def-cbright


# Script to run a SLURM job on a DRAC cluster
#
# This script goes merges a set of sorted files together


filename=$1
output_dir=$2
n=$3

# sorting the files
start=`date +%s`
sort -m $output_dir/*.sorted > $filename.sorted
end=`date +%s`
if [[ $? -eq 0 ]]; then
    echo -e "Merging the file \"$(basename "$filename")\" took $((end - start)) seconds. \n\n" >> "./results/pairs/wts/find_$n/result.log"
fi
