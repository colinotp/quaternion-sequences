#! /bin/bash

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This driver sorts the .pair files generated for a given length n in the first part of the algorithm, using a combination of merge sort and bash quick sort for handling of large files."
	echo "Example usage: ./sortpairs_split.sh n"
	exit 0
fi

if [ ! -e tmp/ ]
then
	mkdir tmp/
fi

n=$1


export LC_ALL=C

# dirname is rowsum_x_y_z_w
for dirname in results/pairs/qts/find_$n/*;
do
	if [ -d $dirname ]
	then
		# filename is pair_AB.pair
		for filename in $dirname/*.pair;
		do
			# Create directories
			echo $filename
			base_name=$(basename "$filename" .pair)

			# output_dir is rowsum_x_y_z_w/pair_AB/
			output_dir="$dirname/$base_name"
			mkdir -p "$output_dir"

			# Split up files into smaller files and sort
			jobids=()
			split -d -l 1000000 "$filename" "$output_dir/${base_name}_part_"

			# file is pair_AB/pair_AB_part_xx
			for file in $output_dir/${base_name}_part_*;
			do
				# Use this for sorting on compute canada
				jobid=$(sbatch job_sort_specific.sh $file $n | awk '{print $4}')
				if [[ -z "$jobid" ]]; then
					echo "Failed to submit job for rowsum $rowsum"
					exit 1
				fi
				jobids+=($jobid)
				echo "Submitted job $jobid for sorting"
				# Use this command for any sorting done not using SLURM
				# sort -S 1G -T tmp/ $filename > $filename.sorted
			done
			# Merge sorted files together
			dep_string=$(IFS=:; echo "${jobids[*]}")
			sbatch --dependency=afterok:$dep_string job_merge_specific.sh $filename $output_dir $n
		done

	fi
done