if [ ! -e tmp/ ]
then
	mkdir tmp/
fi
if [ -z "$1" ]
then
	echo "Need order of PAF files to sort on and optionally case # to solve and compression factor"
	echo "Passing -1 as the case # will solve all cases"
	exit
fi
export LC_ALL=C
n=$1
if [ -n "$3" ]
then
	d=$3
else
	d=1
fi
l=$((n/d))
echo "ORDER $n: Sort compression pairs"
for c in `seq 0 8`
do
	if [ -n "$2" ] && [ $2 -ne $c ] && [ $2 -ne -1 ]
	then
		continue
	fi
	if [ -e matchings/$n.$c.$l.AB.pafs.txt ]
	then
		START=$(date +%s.%N)
		sort -T tmp/ matchings/$n.$c.$l.AB.pafs.txt > matchings/$n.$c.$l.AB.pafs.sorted.txt
		END=$(date +%s.%N)
		DIFF=$(echo "$END - $START" | bc)
		printf "  Case $c: AB PAFs sorted in %.2f seconds\n" $DIFF
		echo $DIFF > timings/$n.$c.$l.AB.sorttime
	fi
	if [ -e matchings/$n.$c.$l.CD.pafs.txt ]
	then
		START=$(date +%s.%N)
		sort -T tmp/ matchings/$n.$c.$l.CD.pafs.txt > matchings/$n.$c.$l.CD.pafs.sorted.txt
		END=$(date +%s.%N)
		DIFF=$(echo "$END - $START" | bc)
		printf "  Case $c: CD PAFs sorted in %.2f seconds\n" $DIFF
		echo $DIFF > timings/$n.$c.$l.CD.sorttime
	fi
done
