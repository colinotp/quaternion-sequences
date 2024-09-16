#! /bin/bash
#SBATCH --time=20:00:00
#SBATCH --mem=8000M
#SBATCH --account=def-cbright
#SBATCH --mail-user=bennet43@uwindsor.ca
#SBATCH --mail-type=BEGIN,FAIL,TIME_LIMIT,TIME_LIMIT_50,TIME_LIMIT_80


filename=$1

# sorting the files
start=`date +%s`
sort -S 1G -T $SLURM_TMPDIR $filename > $filename.sorted
end=`date +%s`
if [[ $? -eq 0 ]]; then
    echo -e "Sorting the file \"$(basename "$filename")\" took $((end - start)) seconds. \n\n" >> "./results/pairs/wts/find_18/result.log"
fi
