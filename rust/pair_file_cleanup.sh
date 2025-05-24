#!/bin/bash

# Cleans up temporary files from given find_n directory (.pair, .pair.sorted, etc.)

n=$1
base_dir="results/pairs/wts"
dir="$base_dir/find_$n"
if [ -d "$dir" ]; then
    rm -rf "$dir"/rowsum_*
fi
