#!/bin/bash
#SBATCH -J FDR
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

Rscript /data/taoyuhuan/projects/FDR/scripts/FDR.R -i $1 -o $2
