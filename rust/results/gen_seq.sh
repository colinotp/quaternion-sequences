#! /bin/sh

# This script generates sequences.
# It is used as follow:
# 
# gen_seq.sh
# to simply run main
# 
# gen_seq.sh type
# to generate all the sequences of a certain type of length n, n in [0,100]
# 
# gen_seq.sh type number
# to generate all the sequences of a certain type of length number
#
#
# Example:
# gen_seq.sh qs 7
# will generate all the perfect quaternion sequences of length 7
#
#
# There are different types:
# qs for perfect quaternion sequences
# ws for Williamson sequences
# wts for Williamson-type sequences


if [ $# -eq 0 ]
then
    cargo run
elif [ $# -eq 1 ]
then
    for i in {1..100}
    do
        ./results/gen_seq.sh $1 $i
    done
else
    filename="results/$1/$2.seq"

    if (test -f "$filename")
    then
        echo "already generated"
        exit 1
    fi

    cargo run $1 $2 > $filename
fi
