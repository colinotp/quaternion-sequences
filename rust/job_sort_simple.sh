#! /bin/bash
#SBATCH --time=100:00:00
#SBATCH --mem=10000M
#SBATCH --account=def-cbright

# "./results/pairs/wts/find_18/rowsum_6_6_0_0/pair_XY.pair"
# "./results/pairs/wts/find_18/rowsum_8_2_2_0/pair_ZX.pair"

sort -S 2G -T $SLURM_TMPDIR "./results/pairs/wts/find_18/rowsum_8_2_2_0/pair_ZX.pair"