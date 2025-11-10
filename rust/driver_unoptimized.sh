#! /bin/bash

# This driver computes sequences for the given length n:
# ./driver.sh qts n
#
# Optional flags:
# -d: delete existing .seq, .pair and .sorted files
# -p <pairing>: run code with chosen pairing of sequences. Options include WX, WY and WZ


if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This driver computes sequences for the given length n, with additional error checking and debugging turned on:"
	echo "./driver_unoptimized.sh <sequencetype> <n> [flags]"
	echo "Optional flags:"
	echo "  * -d: Delete existing .seq, .pair and .sorted files"
	echo "  * -c: Use auto/cross correlation for matching instead of PSD/CPSD"
	echo "  * -p <pairing>: Specify rowsum pairing to be used. Options include XY, XZ, XW. Default is XW"
	echo "  * -h: Convert sequences to Hadamard matrices when finished"
	exit 0
fi

type=$1
n=$2

if [ -z "$type" ] || [ -z "$n" ]; then
	echo 'Incorrect args passed. Try running with --help.'
	exit 1
fi

shift
shift
foldername="./results/pairs/$type/find_$n"
rowsum_pairing="WZ"
hadamard=false
match_option="psd"

# Empty out existing .pair files to avoid conflicts
while getopts "chdp:" flag; do
	case $flag in
		h)
		hadamard=true
		;;
		d)
		./pair_file_cleanup.sh $type $n
		;;
		p)
		rowsum_pairing=$OPTARG
		;;
		c)
		match_option="correlation"
		;;
		/?)
		echo "Invalid argument(s) passed. Exiting."
		exit 1
		;;
	esac
done

# Check if rowsum directories still exist
for d in "$foldername"/rowsum_*; do
  if [ -d "$d" ]; then
    echo "WARNING: results have already been generated for length $n. To run anyway, use the -d flag to overwrite. Exiting."
	exit 1
  fi
done

filename="$foldername/result.log"

if [ ! -e $foldername ]
then
	mkdir -p $foldername
fi

start=`date +%s.%N`

# Creating every necessary file
cargo run pairs $type $n $match_option $rowsum_pairing | tee $filename
if [ $? -ne 0 ]
then
	echo 'ERROR: pairs exited unsuccessfully. See log for additional details'
	exit 1
fi

# Sort pairs
./sortpairs.sh $type $n
if [ $? -ne 0 ]
then
	echo 'ERROR: sort exited unsuccessfully. See log for additional details'
	exit 1
fi

# Matching the file AND reducing to equivalence
cargo run join $type $n | tee $filename -a
if [ $? -ne 0 ]
then
	echo 'ERROR: join exited unsuccessfully. See log for additional details'
	exit 1
fi

if [ $hadamard = true ]; then
	start2=`date +%s.%N`
	cargo run convert hm $type $n | tee $filename -a
	end2=`date +%s.%N`
	elapsed=$(echo "$end2 - $start2" | bc)
	printf "Converting to matrices up to Hadamard equivalence took %.2f seconds\n" $elapsed | tee $filename -a
	
	filename2="$foldername/result.mat"
	matcount=$(wc -l < $filename2)
	echo "$matcount matrices were found after converting up to Hadamard equivalence." | tee $filename -a
fi

echo "Converting PQS to QHM ..." | tee $filename -a
./convert_qhm.sh $type $n

end=`date +%s.%N`
elapsed=$(echo "$end - $start" | bc)
printf "\nTotal execution time was %.2f seconds.\n" $elapsed | tee $filename -a

echo "This output can also be found in $filename"
