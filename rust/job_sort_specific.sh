#! /bin/bash
#SBATCH --time=20:00:00
#SBATCH --mem=8000M
#SBATCH --account=def-cbright

if [ $# -eq 1 ]
then
    filename=$1
else
    filename="./results/pairs/wts/find_17/rowsum_5_5_3_-3/pair_YX.pair"
fi


# sorting the files
start2=`date +%s`
sort -S 2G -T $SLURM_TMPDIR $filename > $filename.sorted
end2=`date +%s`
echo Sorting the files took `expr $end2 - $start2` seconds.
