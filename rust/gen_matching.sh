#! /bin/sh

# This script generates a list of sequences in a .seq file, and a .log file for debug
# It is used as follow: First place yourself in the "rust" folder, then type one of the following commands
# 
# gen_seq.sh
# to simply run main
# 
# gen_seq.sh type
# to generate all the sequences of a certain type of length n, n in [0,100]
# 
# gen_seq.sh type number
# to generate all the sequences of a certain type of length number
#
#
# Example:
# gen_seq.sh qs 7
# will generate all the perfect quaternion sequences of length 7
#
#
# There are different types:
# qs for perfect quaternion sequences
# ws for Williamson sequences
# wts for Williamson-type sequences


number=$1

filename="results/sequences/matches/$number"

cargo run matching $number 2> "$filename.log"
echo "Succesfully generated sequence of type $type and length $number"
exit 1