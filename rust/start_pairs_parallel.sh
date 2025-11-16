#! /bin/bash

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This driver runs the first part of the algorithm with parallelization, stopping after generating the lists of auto/cross correlation values. It is intended for use in an environment with the SLURM job manager."
	echo "Example usage: ./driver_parallel.sh <sequencetype> <n>"
	echo "Optional flags:"
    echo "  * -d: Delete existing .seq, .pair and .sorted files"
	echo "  * -p <pairing>: Specify rowsum pairing to be used. Options include XY, XZ, XW. Default is XW"
	exit 0
fi

type=$1
n=$2
shift
shift
foldername="./results/pairs/$type/find_$n"
rowsum_pairing="XW"
match_option="psd"

while getopts "dcp:" flag; do
	case $flag in
		d)
		./pair_file_cleanup.sh $type $n
		;;
		p)
		rowsum_pairing=$OPTARG
		;;
		c)
		match_option="correlation"
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
./target/release/rust rowsums $type $n 


# read the rowsums file and submit jobs
input="results/pairs/$type/find_$n/rowsums.quad"
while IFS= read -r rowsum
do
    # Submit job for first pair, capturing job ID
	sbatch ./job_pair_single_rowsum.sh $type $n $rowsum $match_option $rowsum_pairing 1
	sbatch ./job_pair_single_rowsum.sh $type $n $rowsum $match_option $rowsum_pairing 2
done < "$input"
