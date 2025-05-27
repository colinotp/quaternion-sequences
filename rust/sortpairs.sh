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


export LC_ALL=C

for dirname in results/pairs/$type/find_$n/*;
do
	if [ -d $dirname ]
	then
		for filename in $dirname/*.pair;
		do
			echo $filename
#			Use this command for sorting on compute canada
#			sort -S 1G -T $SLURM_TMPDIR $filename > $filename.sorted
#			Use this command for any sorting done not using SLURM
			sort -S 1G -T tmp/ $filename > $filename.sorted
		done
	fi
done