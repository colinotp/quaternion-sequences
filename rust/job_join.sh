#! /bin/bash
#SBATCH --time=15:00:00
#SBATCH --mem=4000M
#SBATCH --account=def-cbright


n=17

foldername="./results/pairs/wts/find_$n"
filename="$foldername/result.log"

# sorting the files
start2=`date +%s`
./target/release/rust join $n &>> $filename
end2=`date +%s`
echo Sorting the files took `expr $end2 - $start2` seconds.
