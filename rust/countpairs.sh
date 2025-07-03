#! /bin/bash

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This script counts the number of pairs that were generated for a search of length n:"
	echo "./countpairs.sh sequencetype n"
	exit 0
fi

if [ ! -e tmp/ ]
then
	mkdir tmp/
fi

type=$1
n=$2


count=$((0))

for dirname in results/pairs/$type/find_$n/*;
do
	if [ -d $dirname ]
	then
		for filename in $dirname/*.sorted;
		do
			str=$(wc -l $filename)
			arr=($str)
			nb=${arr[0]}
			count=$(($count + $nb))
		done
	fi
done

echo $count