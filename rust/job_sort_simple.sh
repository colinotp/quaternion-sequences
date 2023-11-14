#! /bin/bash
#SBATCH --time=100:00:00
#SBATCH --mem=10000M
#SBATCH --account=def-cbright


sort -S 2G -T $SLURM_TMPDIR "./results/pairs/wts/find_19/rowsum_3_3_3_7/pair_WY.pair"