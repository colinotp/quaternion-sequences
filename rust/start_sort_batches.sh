#! /bin/bash


n=$1


export LC_ALL=C

for dirname in results/pairs/wts/find_$n/*;
do
	if [ -d $dirname ]
	then
		for filename in $dirname/*.pair;
		do
			echo $filename
			sbatch ./job_sort_specific.sh $filename --output="$filename.sorted"
		done
	fi
done
