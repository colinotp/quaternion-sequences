#! /bin/bash

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This driver exhaustively computes all quaternion-type sequences of length n with parallelization. It is intended for use in an environment with the SLURM job manager."
	echo "Example usage: ./driver.sh n"
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


jobids=()

# read the rowsums file and submit jobs
input="results/pairs/wts/find_$n/rowsums.quad"
while IFS= read -r rowsum
do
	./target/release/rust create $n $rowsum $rowsum_pairing
    # Submit job for first pair, capturing job ID
	jobid=$(sbatch ./job_pair_single_rowsum.sh $n $rowsum $rowsum_pairing 1 | awk '{print $4}')
	# Check for successful job submission
	if [[ -z "$jobid" ]]; then
		echo "Failed to submit job for rowsum $rowsum"
		exit 1
	fi
	jobids+=($jobid)
	echo "Submitted job $jobid for rowsum $rowsum"

	# Submit job for second pair, capturing job ID
	jobid=$(sbatch ./job_pair_single_rowsum.sh $n $rowsum $rowsum_pairing 2 | awk '{print $4}')
	# Check for successful job submission
	if [[ -z "$jobid" ]]; then
		echo "Failed to submit job for rowsum $rowsum"
		exit 1
	fi
	jobids+=($jobid)
	echo "Submitted job $jobid for rowsum $rowsum"

done < "$input"

# Start sorting

# Dependencies for sort jobs
dep_string=$(IFS=:; echo "${jobids[*]}")
jobids2=()

for dirname in results/pairs/wts/find_$n/*;
do
	if [ -d $dirname ]
	then
		for filename in $dirname/*.pair;
		do
			echo $filename
			jobid=$(sbatch --dependency=afterok:$dep_string ./job_sort_specific.sh $filename $n --output="$filename.sorted" | awk '{print $4}')
			# Check for successful job submission
			if [[ -z "$jobid" ]]; then
				echo "Failed to submit sorting job for file $filename"
				exit 1
			fi
			jobids2+=($jobid)
			echo "Submitted job $jobid to sort file $filename"
		done
	fi
done

dep_string2=$(IFS=:; echo "${jobids2[*]}")
sbatch --dependency=afterok:$dep_string2 ./job_join.sh $n

