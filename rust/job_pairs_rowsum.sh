#! /bin/bash
#SBATCH --time=50:00:00
#SBATCH --mem=8000M
#SBATCH --account=def-cbright
#SBATCH --mail-user=bennet43@uwindsor.ca
#SBATCH --mail-type=BEGIN,FAIL,TIME_LIMIT,TIME_LIMIT_50,TIME_LIMIT_80,END

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

foldername="./results/pairs/wts/find_$n"
filename="$foldername/result.log"

if [ ! -e $foldername ]
then
	mkdir $foldername
fi

# Creating every necessary file
start2=`date +%s`
./target/release/rust pairs_rowsum $n $a $b $c $d >> $filename
end2=`date +%s`
echo Creating the sequences took `expr $end2 - $start2` seconds. >> $filename
