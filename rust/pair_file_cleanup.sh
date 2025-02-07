#!/bin/bash

# Cleans up temporary files (.pair, .pair.sorted, etc.)

base_dir="results/pairs/wts"
for i in {1..17}; do
    dir="$base_dir/find_$i"
    if [ -d "$dir" ]; then
        rm -rf "$dir"/rowsum_*
    fi
done