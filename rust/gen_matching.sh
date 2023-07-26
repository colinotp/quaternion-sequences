#! /bin/sh



number=$1

filename="results/sequences/matches/$number"

cargo run matching $number 2> "$filename.log"
echo "Succesfully generated sequence of type $type and length $number"
exit 1