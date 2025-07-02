#! /bin/bash

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This script submits a job to SLURM for each .pair file corresponding to a given length n, and sorts them."
	echo "Example usage: ./start_sort_batches <sequencetype> <n>"
	exit 0
fi

type=$1
n=$2


export LC_ALL=C

for dirname in results/pairs/$type/find_$n/*;
do
	if [ -d $dirname ]
	then
		for filename in $dirname/*.pair;
		do
			echo $filename
			sbatch ./job_sort_specific.sh $type $filename $n --output="$filename.sorted"
		done
	fi
done
