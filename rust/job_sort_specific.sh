#! /bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=4000M
#SBATCH --account=def-cbright


filename=$1
n=$2

# sorting the files
start=`date +%s`
sort -S 2G -T $SLURM_TMPDIR $filename > $filename.sorted
end=`date +%s`
if [[ $? -eq 0 ]]; then
    echo -e "Sorting the file \"$(basename "$filename")\" took $((end - start)) seconds. \n\n" >> "./results/pairs/wts/find_$n/result.log"
fi
