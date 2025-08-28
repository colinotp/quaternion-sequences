#! /bin/bash

# This driver computes sequences for the given length n:
# ./driver.sh qts n
#
# Optional flags:
# -d: delete existing .seq, .pair and .sorted files
# -p <pairing>: run code with chosen pairing of sequences. Options include XY, XW, XZ


if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This driver computes sequences for the given length n. It does not convert the sequences to Hadamard matrices unless the -h flag is passed. Usage:"
	echo "./driver.sh <sequencetype> <n> [flags]"
	echo "Optional flags:"
	echo "  * -s: Use this flag for SLURM jobs"
	echo "  * -h: Convert sequences to Hadamard matrices when finished"
	echo "  * -d: Delete existing .seq, .pair and .sorted files"
	echo "  * -c: Use auto/cross correlation for matching instead of PSD/CPSD"
	echo "  * -p <pairing>: Specify rowsum pairing to be used. Options include WX, WY and WZ (e.g., WX means that the sequences of rowsum W are paired with the sequences of rowsum X). Note that the code follows the convention W <= X <= Y <= Z. Default is WZ"
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
use_slurm=false
hadamard=false
match_option="psd"

# Empty out existing .pair files to avoid conflicts
while getopts "chsdp:" flag; do
	case $flag in
		s)
		use_slurm=true
		;;
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

release_dir="./target/release/rust"
if [ $use_slurm = false ]; then
	# Cargo automatically checks if source files have changed before deciding whether to compile
	cargo build --release
elif [ ! -e $release_dir ]; then
	echo "ERROR: Binary not found, compile with 'cargo build --release'"
fi

filename="$foldername/result.log"

if [ ! -e $foldername ]
then
	mkdir -p $foldername
fi

start=`date +%s.%N`

# Creating every necessary file
./target/release/rust pairs $type $n $match_option $rowsum_pairing | tee $filename
if [ $? -ne 0 ]
then
	echo 'ERROR: pairs exited unsuccessfully. See log for additional details'
	exit 1
fi

# sorting the files
start2=`date +%s.%N`
if [ "$use_slurm" = true ]; then
	./sortpairs.sh $type $n -s | tee $filename -a
else
	./sortpairs.sh $type $n | tee $filename -a
fi

if [ $? -ne 0 ]
then
	echo 'ERROR: sort exited unsuccessfully. See log for additional details'
	exit 1
fi
end2=`date +%s.%N`
elapsed=$(echo "$end2 - $start2" | bc)
printf "Sorting the files took: %.2f seconds\n\n" $elapsed
printf "Sorting the files took: %.2f seconds\n\n" $elapsed >> $filename

# Matching the file AND reducing to equivalence
start2=`date +%s`
./target/release/rust join $type $n | tee $filename -a
if [ $? -ne 0 ]
then
	echo 'ERROR: join exited unsuccessfully. See log for additional details'
	exit 1
fi
end2=`date +%s`


if [ $hadamard = true ]; then
	start2=`date +%s.%N`
	./target/release/rust convert $type $n | tee $filename -a
	end2=`date +%s.%N`start2=`date +%s.%N`
	elapsed=$(echo "$end2 - $start2" | bc)
	printf "Converting to matrices up to Hadamard equivalence took %.2f seconds\n" $elapsed
	printf "Converting to matrices up to Hadamard equivalence took %.2f seconds\n" $elapsed >> $filename
	
	filename2="$foldername/result.mat"
	matcount=$(wc -l < $filename2)
	echo "$matcount matrices were found after converting up to Hadamard equivalence." | tee $filename -a
fi

end=`date +%s.%N`
elapsed=$(echo "$end - $start" | bc)
printf "Total execution time was %.2f seconds.\n" $elapsed
printf "Total execution time was %.2f seconds.\n" $elapsed >> $filename

echo "This output can also be found in $filename"
