#! /bin/sh


if [ $# -eq 0 ] || [ "$1" = "help" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]
then
	echo "This driver computes the perfect Williamson Type sequences for the given n:"
	echo "./driver.sh n"
	exit 0
fi

n=$1

start=`date +%s`

cargo run pairs $n
./sortpairs.sh wts $n
cargo run join $n

end=`date +%s`
echo Execution time was `expr $end - $start` seconds.