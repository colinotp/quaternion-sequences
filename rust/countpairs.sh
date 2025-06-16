#! /bin/bash

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This script counts the number of pairs that were generated for a search of length n:"
	echo "./countpairs.sh n"
	exit 0
fi

if [ ! -e tmp/ ]
then
	mkdir tmp/
fi

n=$1


count=$((0))

for dirname in results/pairs/qts/find_$n/*;
do
	if [ -d $dirname ]
	then
		for filename in $dirname/*.pair;
		do
			str=$(wc -l $filename)
			arr=($str)
			nb=${arr[0]}
			count=$(($count + $nb))
		done
	fi
done

echo $count