

n=$1

for filename in "result/pairs/*.pair"
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
