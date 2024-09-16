#! /bin/sh

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
			sort -S 1G -T $SLURM_TMPDIR $filename > $filename.sorted
		done
	fi
done