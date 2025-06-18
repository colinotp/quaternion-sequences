#! /bin/bash
#SBATCH --time=300:00:00
#SBATCH --mem=70000M
#SBATCH --account=def-cbright
#SBATCH --cpus-per-task=16

# Script to run a SLURM job on a DRAC cluster
#
# This script converts found sequences into Hadamard matrices (up to Hadamard equivalence)

type=$1
n=$2


# sorting the files
start2=`date +%s`
./target/release/rust convert $type $n
end2=`date +%s`
echo Converting to matrices up to Hadamard equivalence took `expr $end2 - $start2` seconds.
