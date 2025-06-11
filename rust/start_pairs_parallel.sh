#! /bin/bash

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This driver runs the first part of the algorithm with parallelization, stopping after generating the lists of auto/cross correlation values. It is intended for use in an environment with the SLURM job manager."
	echo "Example usage: ./driver_parallel.sh n"
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

# Generate rowsums
./target/release/rust rowsums $n 


# read the rowsums file and submit jobs
input="results/pairs/wts/find_$n/rowsums.quad"
while IFS= read -r rowsum
do
    # Submit job for first pair, capturing job ID
	sbatch ./job_pair_single_rowsum.sh $n $rowsum $rowsum_pairing 1
	sbatch ./job_pair_single_rowsum.sh $n $rowsum $rowsum_pairing 2
done < "$input"
