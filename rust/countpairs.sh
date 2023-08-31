#! /bin/sh

if [ ! -e tmp/ ]
then
	mkdir tmp/
fi

n=$1


count=$((0))

for dirname in results/pairs/wts/find_$n/*;
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