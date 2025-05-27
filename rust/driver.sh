#! /bin/bash

# This driver computes the perfect Williamson Type sequences for the given n:
# ./driver.sh n
#
# At the moment, it doesn't return the list of sequences, but just the number of sequences.
# I'll add it quickly.


if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This driver computes the perfect Williamson Type sequences for the given n:"
	echo "./driver.sh n"
	exit 0
fi


n=$1
shift
foldername="./results/pairs/wts/find_$n"
rowsum_pairing="XW"

# Empty out existing .pair files to avoid conflicts
while getopts "dp:" flag; do
	case $flag in
		d)
		./pair_file_cleanup.sh $n
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
	mkdir $foldername
fi

start=`date +%s`

# Creating every necessary file
start2=`date +%s`
./target/release/rust pairs $n $rowsum_pairing &> $filename
end2=`date +%s`
echo Creating the sequences took `expr $end2 - $start2` seconds. 
echo -e Creating the sequences took `expr $end2 - $start2` seconds. "\n \n" >> $filename

# sorting the files
start2=`date +%s`
./sortpairs.sh wts $n &>> $filename
end2=`date +%s`
echo Sorting the files took `expr $end2 - $start2` seconds.
echo -e Sorting the files took `expr $end2 - $start2` seconds. "\n \n" >> $filename

# Matching the file AND reducing to equivalence
start2=`date +%s`
./target/release/rust join $n &>> $filename
end2=`date +%s`
echo Matching the sequences took `expr $end2 - $start2` seconds.
echo -e Matching the sequences took `expr $end2 - $start2` seconds. "\n \n" >> $filename

end=`date +%s`
echo Total execution time was `expr $end - $start` seconds.
echo -e Total execution time was `expr $end - $start` seconds. >> $filename