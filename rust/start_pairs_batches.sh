#! /bin/bash

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This script submits a job to SLURM for each rowsum corresponding to a given length n."
    echo "The job executes the first part of the algorithm, and stops after generating the lists of auto and cross correlation values for the pairs."
	echo "Example usage: ./start_pairs_batches n -p XY"
    echo "The -p flag is an optional flag that allows you to specify the pairing of the sequences based off of their rowsum."
	exit 0
fi

n=$1
shift
rowsum_pairing="XW"

while getopts "dp:" flag; do
	case $flag in
		d)
		./pair_file_cleanup.sh $n
		;;
		p)
		rowsum_pairing=$OPTARG
		;;
		/?)
		echo "Invalid argument(s) passed. Exiting."
		exit 1
		;;
	esac
done

# Check if rowsum directories still exist
for d in "$foldername"/rowsum_*; do
  if [ -d "$d" ]; then
    echo "WARNING: results have already been generated for length $n. To run anyway, use the -d flag to overwrite. Exiting."
	exit 1
  fi
done

./target/release/rust rowsums $n 



# read the rowsums file
input="results/pairs/wts/find_$n/rowsums.quad"
while IFS= read -r rowsum
do
    #launch the batches for each rowsum
    sbatch ./job_pairs_rowsum.sh $n $rowsum $rowsum_pairing
done < "$input"
