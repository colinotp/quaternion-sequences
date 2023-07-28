#! /bin/sh

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

start=`date +%s`

# Creating every necessary file
start2=`date +%s`
cargo run pairs $n
end2=`date +%s`
echo Sorting the files took `expr $end2 - $start2` seconds.

# sorting the files
start2=`date +%s`
./sortpairs.sh wts $n
end2=`date +%s`
echo Sorting the files took `expr $end2 - $start2` seconds.

# Matching the file AND reducing to equivalence
start2=`date +%s`
cargo run join $n
end2=`date +%s`
echo Sorting the files took `expr $end2 - $start2` seconds.

end=`date +%s`
echo Execution time was `expr $end - $start` seconds.