#! /bin/sh

# WIP
# Intended to split up the files being sorted into smaller files, and sort them individually to ease filesystem load

if [ ! -e tmp/ ]
then
	mkdir tmp/
fi

type=$1
n=$2


export LC_ALL=C

for dirname in results/pairs/$type/find_$n/*;
do
	if [ -d $dirname ]
	then
		for filename in $dirname/*.pair;
		do
			echo $filename
			for filename in "$dirname"/*.pair;
			do
				base_name=$(base_name "$filename" .pair)

				output_dir="$dirname/$base_name"
				mkdir -p "$output_dir"

				split -d -l 15000000 "$filename" "$output_dir/${base_name}_part_"
			done
		done
	fi
done