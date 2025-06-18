#!/bin/bash

# Cleans up temporary files from find_n directories (.pair, .pair.sorted, etc.)

type=$1

if [ $# -eq 2 ]; then
    START=$2
    END=$2
elif [ $# -eq 3 ]; then
    START=$2
    END=$3
else
    echo "This script deletes rowsum_w_x_y_z directories (temporary files containing pair generation). Usage:"
    echo "  * To clean directories corresponding to sequences of type 'seqtype' for a single length n, use ./pair_file_cleanup.sh seqtype n"
    echo "  * To clean directories corresponding to sequences of type 'seqtype' for all lengths from a to b (inclusive), use ./pair_file_cleanup.sh seqtype a b"
    exit 0
fi

base_dir="results/pairs/$type"
for (( i=START; i<=END; i++ ))
do
    dir="$base_dir/find_$i"
    if [ -d "$dir" ]; then
        rm -rf "$dir"/rowsum_*
    fi
done
