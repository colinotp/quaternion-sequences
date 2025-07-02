#! /bin/bash

# This driver computes sequences for the given length n:
# ./driver_unoptimized.sh qts n
#
# Optional flags:
# -d: delete existing .seq, .pair and .sorted files
# -p <pairing>: run code with chosen pairing of sequences. Options include XY, XW, XZ


if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This driver computes sequences for the given length n, with additional error checking and debugging turned on:"
	echo "./driver.sh <sequencetype> <n>"
	echo "Optional flags:"
	echo "  * -d: Delete existing .seq, .pair and .sorted files"
	echo "  * -p <pairing>: Specify rowsum pairing to be used. Options include XY, XZ, XW. Default is XW"
	exit 0
fi

type=$1
n=$2
shift
shift
foldername="./results/pairs/$type/find_$n"
rowsum_pairing="XW"

# Empty out existing .pair files to avoid conflicts
while getopts "dp:" flag; do
	case $flag in
		d)
		./pair_file_cleanup.sh $type $n
		;;
		p)
		rowsum_pairing=$OPTARG
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

start=`date +%s`

# Creating every necessary file
start2=`date +%s`
cargo run pairs $type $n $rowsum_pairing &> $filename
if [ $? -ne 0 ]
then
	echo 'ERROR: pairs exited unsuccessfully. See log for details'
	exit 1
fi
end2=`date +%s`
echo Creating the sequences took `expr $end2 - $start2` seconds. 
echo -e Creating the sequences took `expr $end2 - $start2` seconds. "\n \n" >> $filename

# sorting the files
start2=`date +%s`
./sortpairs.sh $type $n &>> $filename
if [ $? -ne 0 ]
then
	echo 'ERROR: sort exited unsuccessfully. See log for details'
	exit 1
fi
end2=`date +%s`
echo Sorting the files took `expr $end2 - $start2` seconds.
echo -e Sorting the files took `expr $end2 - $start2` seconds. "\n \n" >> $filename

# Matching the file AND reducing to equivalence
start2=`date +%s`
cargo run join $type $n &>> $filename
if [ $? -ne 0 ]
then
	echo 'ERROR: join exited unsuccessfully. See log for details'
	exit 1
fi
end2=`date +%s`
echo Matching the sequences took `expr $end2 - $start2` seconds.
echo -e Matching the sequences took `expr $end2 - $start2` seconds. "\n \n" >> $filename

end=`date +%s`
echo Total execution time was `expr $end - $start` seconds.
echo -e Total execution time was `expr $end - $start` seconds. >> $filename
