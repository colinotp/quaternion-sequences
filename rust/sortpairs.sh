#! /bin/sh

type=$1
n=$2

for dirname in results/pairs/$type/find_$n/*;
do
	for filename in $dirname/*.pair;
	do
		echo $filename
		sort -g $filename > $filename.sorted
	done
done
