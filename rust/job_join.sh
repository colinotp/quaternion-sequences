#! /bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=4000M
#SBATCH --account=def-cbright

# Script to run a SLURM job on a DRAC cluster
#
# This script runs the part of the algorithm that goes
# through the sorted auto and cross correlation values
# to find valid QTS, and then computes the corresponding PQS

type=$1
n=$2

foldername="./results/pairs/qts/find_$n"
filename="$foldername/result.log"

# sorting the files
start2=`date +%s`
./target/release/rust join $type $n &>> $filename
end2=`date +%s`
if [[ $? -eq 0 ]]; then
    echo -e "Joining the files together took $((end2 - start2)) seconds. \n\n" >> "./results/pairs/$type/find_$n/result.log"
fi

