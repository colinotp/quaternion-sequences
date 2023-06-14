#!/bin/bash

if [ -z "$1" ]
then
	echo "Use ./driver.sh n to enumerate quasi-Williamson matrices of orders n."
	echo "Optionally, pass the compression factor d as the second argument (d must divide n).  By default no compression is used."
	echo "Warning: Currently does not consider all rowsum decompositions and therefore the search will not be exhaustive."
	exit
else
	make all
	n=$1
	if [ -n "$2" ]
	then
		d=$2
	else
		d=1
	fi
	./generate_matching_instances $n
	./generate_compression_matchings $n $d
	./generate_pairedmatchings $n $d
	./sortpairs.sh $n -1 $d
	./join_pairedmatchings $n -1 $d
	./remove_equivalent_matchedpairs $n -1 $d
	if (( d > 1 ))
	then
		if [ ! -f maplesat_static ]
		then
			./compile_maplesat.sh
		fi
		python2 generate_compstring_instances.py $n
		python2 solve_compstring_instances.py $n -1 $d
	else
		mkdir -p exhaust
		c=0
		while [ -f matchedseqns/$n.$c.$n.inequiv ]
		do
			cp matchedseqns/$n.$c.$n.inequiv exhaust/$n.$c
			c=$((c+1))
		done
	fi
	./remove_equivalent_exhaust $n
fi
