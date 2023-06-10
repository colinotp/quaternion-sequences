#! /bin/sh

# This script generates a list of sequences in a .seq file, and a .log file for debug
# It is used as follow: First place yourself in the "rust" folder, then type one of the following commands
# 
# results/gen_seq.sh
# to simply run main
# 
# results/gen_seq.sh type
# to generate all the sequences of a certain type of length n, n in [0,100]
# 
# results/gen_seq.sh type number
# to generate all the sequences of a certain type of length number
#
#
# Example:
# results/gen_seq.sh qs 7
# will generate all the perfect quaternion sequences of length 7
#
#
# There are different types:
# qs for perfect quaternion sequences
# ws for Williamson sequences
# wts for Williamson-type sequences


type=$1
number=$2
regen=0

if [ $# -gt 0 ]
then #check for the --regen keyword
    if [ $1 = "--regen" ]
    then
        type=$2
        number=$3
        regen=1
    fi
fi


if [ $# -eq $regen ]
then # run main
    cargo run
elif [ $# -eq $(($regen+1)) ]
then # generate all the sequences of a certain type of length n, n in [0,100]
    for i in {1..100}
    do
        if [ $regen -eq 1 ]
        then
            ./results/gen_seq.sh --regen $type $i
        else
            ./results/gen_seq.sh $type $i
        fi
    done
else # generate all the sequences of a certain type of length number
    filename="results/$type/$number"

    if (test -f "$filename.log")
    then
        if ! (cat "$filename.log" | grep -q -e "error" || [ $regen -eq 1 ] )
        then
            echo "already generated sequence of type $type and length $number"
            exit 1
        fi
    fi

    cargo run $type $number > "$filename.seq" 2> "$filename.log"
    echo "Succesfully generated sequence of type $type and length $number"
    exit 1
    
fi
