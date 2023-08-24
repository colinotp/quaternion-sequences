#! /bin/bash
#SBATCH --time=80:00:00
#SBATCH --mem=4000M
#SBATCH --account=def-cbright

if [ $# -ne 5 ]
then
    echo "not enough arguments"
    exit 1
fi

n=$1
a=$2
b=$3
c=$4
d=$5
start=`date +%s`


# go through rowsums
# start all the batches


# Creating every necessary file
start2=`date +%s`
./target/release/rust pairs_rowsum $n $a $b $c $d
end2=`date +%s`
echo Creating the sequences took `expr $end2 - $start2` seconds. 
