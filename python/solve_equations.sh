#! /bin/sh





type=$1
number=$2


foldername="results/equations/$type/find_$number"
for folder in ../rust/$foldername/*
do
    for file in $folder/*
    do
        if ! (test -f $file)
        then
            echo "$file"
            echo "File doesn't exist"
            exit 1
        fi

        python3 gurobi.py $file

        if (test -f "result.sol")
        then
            echo $file
            exit 0
        fi

    done
done
