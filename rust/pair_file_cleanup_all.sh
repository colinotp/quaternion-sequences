#!/bin/bash

# Cleans up temporary files from all find_n directories (.pair, .pair.sorted, etc.)

base_dir="results/pairs/wts"
for i in {1..19}; do
    ./pair_file_cleanup.sh $i
done
