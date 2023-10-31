#! /bin/bash
#SBATCH --time=150:00:00
#SBATCH --mem=16000M
#SBATCH --account=def-cbright

n=$1


# sorting the files
start2=`date +%s`
./target/release/rust convert $n
end2=`date +%s`
echo Converting to matrices up to Hadamard equivalence took `expr $end2 - $start2` seconds.
