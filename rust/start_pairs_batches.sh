#! /bin/bash


n=$1
./target/release/rust rowsums $n 


# read the rowsums file
input="results/pairs/wts/find_$n/rowsums.quad"
while IFS= read -r rowsum
do
    #launch the batches for each rowsum
    sbatch ./job_pairs_rowsum.sh $n $rowsum
done < "$input"
