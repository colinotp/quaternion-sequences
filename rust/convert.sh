#! /bin/bash
# This script converts found sequences into Hadamard matrices (up to Hadamard equivalence)

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This script converts found sequences into Hadamard matrices (up to Hadamard equivalence). Usage:"
	echo "./convert.sh <sequencetype> <length>"
	exit 0
fi

type=$1
n=$2


# sorting the files
start2=`date +%s`
./target/release/rust convert $type $n
end2=`date +%s`
echo Converting to matrices up to Hadamard equivalence took `expr $end2 - $start2` seconds.
