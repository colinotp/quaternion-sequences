#! /bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --mem=8000M
#SBATCH --account=def-cbright

# Script to run a SLURM job on a DRAC cluster
#
# This script runs the first part of the algorithm for a single set of rowsums
# Stops after generating the lists of auto and cross correlation values for the pairs
# Normally called from start_pairs_parallel.sh

if [ $# -ne 8 ]
then
    echo "not enough arguments"
    exit 1
fi

type=$1
n=$2
a=$3
b=$4
c=$5
d=$6
match_option=$7
rowsum_pairing=$8
pair=$9
start=`date +%s`


# go through rowsums
# start all the batches

foldername="./results/pairs/$type/find_$n"
filename="$foldername/result.log"

if [ ! -e $foldername ]
then
	mkdir $foldername
fi

# Creating every necessary file
start2=`date +%s`
./target/release/rust pair_single $type $n $a $b $c $d $match_option $rowsum_pairing $pair &>> $filename
end2=`date +%s`
echo Creating the sequences took `expr $end2 - $start2` seconds. >> $filename
