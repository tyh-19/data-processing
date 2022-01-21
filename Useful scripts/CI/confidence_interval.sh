#!/bin/bash
#SBATCH -J confidence_interval
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

Rscript /data/taoyuhuan/projects/exOmics_RNA/level_3_Insert_length/Confidence_Interval.R -i /data/taoyuhuan/projects/exOmics_RNA/level_3_Insert_length/$1/summary/miso_Insert_length_summary.txt -o /data/taoyuhuan/projects/exOmics_RNA/level_3_Insert_length/$1/summary/length_CI_median.txt
