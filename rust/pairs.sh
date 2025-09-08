#! /bin/bash

if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This driver runs the first part of the algorithm for a given length n, generating rowsum decompositions and creating .pair files with PSD/CPSD data. Usage:"
	echo "./pairs.sh <sequencetype> <n> [flags]"
	echo "Optional flags:"
    echo "  * -s: Use this flag for SLURM jobs"
	echo "  * -d: Delete existing .seq, .pair and .sorted files"
	echo "  * -c: Use auto/cross correlation for matching instead of PSD/CPSD"
	echo "  * -p <pairing>: Specify rowsum pairing to be used. Options include WX, WY and WZ (e.g., WX means that the sequences of rowsum W are paired with the sequences of rowsum X). Note that the code follows the convention W <= X <= Y <= Z. Default is WZ"
	exit 0
fi

type=$1
n=$2

# Folder to store output
foldername="./results/pairs/$type/find_$n"

if [ -z "$type" ] || [ -z "$n" ]; then
	echo 'Incorrect args passed. Try running with --help.'
	exit 1
fi

# Manage flags
shift
shift

use_slurm=false
rowsum_pairing="WZ"
match_option="psd"
while getopts "cdsp:" flag; do
	case $flag in
        s)
		use_slurm=true
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

# If not part of a slurm job, we can compile without user input
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


# Call rust code
start=`date +%s`
./target/release/rust pairs $type $n $match_option $rowsum_pairing | tee $filename
end=`date +%s`
echo Generating the .pair files took `expr $end - $start` seconds. 
echo -e Generating the .pair files took `expr $end - $start` seconds. "\n \n" >> $filename
