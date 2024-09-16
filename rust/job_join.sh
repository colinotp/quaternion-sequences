#! /bin/bash
#SBATCH --time=30:00:00
#SBATCH --mem=4000M
#SBATCH --account=def-cbright
#SBATCH --mail-user=bennet43@uwindsor.ca
#SBATCH --mail-type=BEGIN,FAIL,TIME_LIMIT,TIME_LIMIT_50,TIME_LIMIT_80,END


n=18

foldername="./results/pairs/wts/find_$n"
filename="$foldername/result.log"

# sorting the files
start2=`date +%s`
./target/release/rust join $n &>> $filename
end2=`date +%s`
if [[ $? -eq 0 ]]; then
    echo -e "Joining the files together took $((end - start)) seconds. \n\n" >> "./results/pairs/wts/find_18/result.log"
fi

