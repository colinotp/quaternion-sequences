#! /bin/bash

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This script submits a job to SLURM for each rowsum corresponding to a given length n."
    echo "The job executes the first part of the algorithm, and stops after generating the lists of auto and cross correlation values for the pairs."
	echo "Example usage: ./start_pairs_batches n"
	exit 0
fi

n=$1
./target/release/rust rowsums $n 


# read the rowsums file
input="results/pairs/wts/find_$n/rowsums.quad"
while IFS= read -r rowsum
do
    #launch the batches for each rowsum
    sbatch ./job_pairs_rowsum.sh $n $rowsum
done < "$input"
