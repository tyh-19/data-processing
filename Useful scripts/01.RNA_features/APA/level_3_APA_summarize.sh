#!/bin/bash
#SBATCH -J APA
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
positive=$2
negative=$3
input="/data/taoyuhuan/projects/exOmics_RNA/${dataset}/output/bam"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_APA/${dataset}"
mkdir -p ${output}/wig
mkdir -p ${output}/dapar
mkdir -p ${output}/config
mkdir -p ${output}/matrix
mkdir -p ${output}/matrix/${positive}vs${negative}

echo "Start summarize ${dataset}(positive:${positive},negative:${negative}) at `date`"
#Summarize
/data/taoyuhuan/projects/exOmics_RNA/bin/scripts/parseDapar.py -i ${output}/dapar/${positive}vs${negative}/result_All_Prediction_Results.txt -c ${output}/config/config_${positive}vs${negative}.txt -l ${output}/matrix/${positive}vs${negative}/long.txt -s ${output}/matrix/${positive}vs${negative}/short.txt -p ${output}/matrix/${positive}vs${negative}/PDUI.txt

echo "Done!"
