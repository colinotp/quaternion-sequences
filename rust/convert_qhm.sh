#! /bin/bash
# This script converts found sequences into Hadamard matrices (up to Hadamard equivalence)

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This script converts found PQS into Quaternionic Hadamard matrices. Usage:"
	echo "./convert_qhm.sh <sequencetype> <length>"
	exit 0
fi

type=$1
n=$2
filename="./results/pairs/$type/find_$n/result.log"

# sorting the files
start2=`date +%s.%N`
./target/release/rust convert qhm $type $n | tee -a $filename
end2=`date +%s.%N`
elapsed=$(echo "$end2 - $start2" | bc)
printf "Converting to PQS to QHM took %.2f seconds\n" $elapsed
printf "Converting to PQS to QHM took %.2f seconds\n" $elapsed >> $filename
