#! /bin/bash
#SBATCH --time=00:05:00
#SBATCH --mem=2000M
#SBATCH --account=def-cbright

if [[ $# -ne 1 ]]; then
    echo 'Must pass the size'
fi


n=$1

./target/release/rust rowsums 18



