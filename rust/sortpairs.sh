#! /bin/sh

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This script sorts all of the generated .pair files for a given length n:"
	echo "./driver.sh <sequencetype> <n>"
	echo "Optional flags:"
	echo "  * -s: Use when sorting in a SLURM job"
	exit 0
fi

if [ ! -e tmp/ ]
then
	mkdir tmp/
fi

type=$1
n=$2
shift
shift

results="./results/pairs/$type/find_$n/result.log"
use_slurm=false
# Empty out existing .pair files to avoid conflicts
while getopts "s" flag; do
	case $flag in
		s)
		use_slurm=true
		;;
		/?)
		echo "Invalid argument(s) passed. Exiting."
		exit 1
		;;
	esac
done

export LC_ALL=C

start=`date +%s.%N`
for dirname in results/pairs/$type/find_$n/*;
do
	if [ -d $dirname ]
	then
		for filename in $dirname/*.pair;
		do
			# For consistency use a single core to sort, even when more cores are available
			echo "Sorting file $filename ..." | tee $results -a
			if [ "$use_slurm" = true ]; then
				sort --parallel=1 -S 1G -T $SLURM_TMPDIR $filename > $filename.sorted
				status=$?
			else
				sort --parallel=1 -S 1G -T tmp/ $filename > $filename.sorted
				status=$?
			fi

			# If sort was successful, remove .pair files
			if [ $status -eq 0 ]; then
				rm $filename
			fi
		done
	fi
done

end=`date +%s.%N`
elapsed=$(echo "$end - $start" | bc)
printf "Total time to sort: %.2f seconds.\n\n" $elapsed | tee $results -a
