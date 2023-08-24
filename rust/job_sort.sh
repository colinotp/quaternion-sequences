#! /bin/bash
#SBATCH --time=3:00:00
#SBATCH --mem=4000M
#SBATCH --account=def-cbright


# sorting the files
start2=`date +%s`
./sortpairs.sh wts 17
end2=`date +%s`
echo Sorting the files took `expr $end2 - $start2` seconds.
