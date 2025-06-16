# !/bin/bash

if [ $# -eq 1 ]; then
    START=$1
    END=$1
elif [ $# -eq 2 ]; then
    START=$1
    END=$2
else
    echo "This script deletes rowsum_w_x_y_z directories (temporary files containing pair generation). Usage:"
    echo "  * To clean directories corresponding to a single length n, use ./pair_file_cleanup.sh n"
    echo "  * To clean directories corresponding to all lengths from a to b (inclusive), use ./pair_file_cleanup.sh a b"
    exit 0
fi

base_dir="results/pairs/qts"
for (( i=START; i<=END; i++ ))
do
    ./target/release/rust amicable $i
done
