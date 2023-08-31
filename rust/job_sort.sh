#! /bin/bash
#SBATCH --time=30:00:00
#SBATCH --mem=4000M
#SBATCH --account=def-cbright

n=$1

# sorting the files
start2=`date +%s`
./sortpairs.sh wts $n
end2=`date +%s`
echo Sorting the files took `expr $end2 - $start2` seconds.
