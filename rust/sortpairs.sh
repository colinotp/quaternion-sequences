#! /bin/sh

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This script sorts all of the generated .pair files for a given length n:"
	echo "./driver.sh wts n"
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

for dirname in results/pairs/$type/find_$n/*;
do
	if [ -d $dirname ]
	then
		for filename in $dirname/*.pair;
		do
			echo $filename
			if [ "$use_slurm" = true ]; then
				sort -S 1G -T $SLURM_TMPDIR $filename > $filename.sorted
			else
			sort -S 1G -T tmp/ $filename > $filename.sorted
			fi
		done
	fi
done