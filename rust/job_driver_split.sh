#!/bin/bash

if [[ $# -ne 1 ]]; then
	echo 'Incorrect arguments. Correct usage:'
	echo './job_driver_split.sh <length>'
	exit 1
fi

n=$1

foldername="./results/pairs/wts/find_$n"
filename="$foldername/result.log"

if [ ! -e $foldername ]
then
	mkdir $foldername
fi

start=`date +%s`

echo 'Computing rowsums...'
./target/release/rust rowsums 18
echo 'Rowsums computed successfully'

rowsumfile="./results/pairs/wts/find_$n/rowsums.quad"

if [[ ! -f "$rowsumfile" ]]; then
    echo "File $rowsumfile not found"
    exit 1
fi

while IFS= read -r line; do
	set -- $line


done < "$rowsumfile"

end=`date +%s`
echo Total execution time was `expr $end - $start` seconds.
echo -e Total execution time was `expr $end - $start` seconds. >> $filename

# Need to call main with arg `rowsums` to generate rowsums
# Then call job_pairs_rowsums.sh with args `length, r1,..,r4`
# Need to figure out how much that accomplishes and what needs to be done after
# Also need to make sure we are ending jobs after each sort if possible to clear $SLURM_TMPDIR