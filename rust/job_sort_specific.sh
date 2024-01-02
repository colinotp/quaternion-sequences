#! /bin/bash
#SBATCH --time=20:00:00
#SBATCH --mem=8000M
#SBATCH --account=def-cbright

# sorting the files
sort -S 2G -T $SLURM_TMPDIR $1
