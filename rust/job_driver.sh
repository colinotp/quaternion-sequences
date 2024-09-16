#! /bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=2000M
#SBATCH --account=def-cbright
#SBATCH --mail-user=bennet43@uwindsor.ca
#SBATCH --mail-type=BEGIN,FAIL,TIME_LIMIT,TIME_LIMIT_50,TIME_LIMIT_80,END

./driver.sh 17
