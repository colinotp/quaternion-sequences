#! /bin/bash


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
    # Submit job for first pair, capturing job ID
	jobid=$(sbatch ./job_pair_single_rowsum.sh $n $rowsum $rowsum_pairing 1 | awk '{print $4}')
	# Check for successful job submission
	if [[ -z "$jobid" ]]; then
		echo "Failed to submit job for rowsum $rowsum"
		exit 1
	fi
	jobids+=($jobid)
	echo "Submitted job $jobid for rowsum $rowsum"

	# Submit job for first pair, capturing job ID
	jobid=$(sbatch ./job_pair_single_rowsum.sh $n $rowsum $rowsum_pairing 2 | awk '{print $4}')
	# Check for successful job submission
	if [[ -z "$jobid" ]]; then
		echo "Failed to submit job for rowsum $rowsum"
		exit 1
	fi
	jobids+=($jobid)
	echo "Submitted job $jobid for rowsum $rowsum"

done < "$input"

dep_string=$(IFS=:; echo "${jobids[*]}")
jobid=$(sbatch --dependency=afterok:$dep_string ./job_sort.sh wts $n | awk '{print $4}')
# Check for successful job submission
if [[ -z "$jobid" ]]; then
	echo "Failed to submit sort job"
	exit 1
fi
echo "Submitted job $jobid for sorting"

sbatch --dependency=afterok:$jobid ./job_join.sh $n

