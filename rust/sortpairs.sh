#! /bin/sh

type=$1
n=$2


export LC_ALL=C

for dirname in results/pairs/$type/find_$n/*;
do
	for filename in $dirname/*.pair;
	do
		echo $filename
		sort $filename > $filename.sorted
	done
done
