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

for dirname in results/pairs/wts/find_$n/*;
do
	if [ -d $dirname ]
	then
		for filename in $dirname/*.pair;
		do
			# Create directories
			echo $filename
			base_name=$(basename "$filename" .pair)
			output_dir="$dirname/$base_name"
			mkdir -p "$output_dir"

			# Split up files into smaller files and sort
			split -d -l 1000000 "$filename" "$output_dir/${base_name}_part_"
			for file in $output_dir/${base_name}_part_*;
			do
				# Use this command for sorting on compute canada
				sort -S 1G -T $SLURM_TMPDIR $filename > $filename.sorted
				# Use this command for any sorting done not using SLURM
				# sort -S 1G -T tmp/ $filename > $filename.sorted
			done
			# Merge sorted files together
			sort -m $output_dir/*.sorted > $filename.sorted
		done

	fi
done