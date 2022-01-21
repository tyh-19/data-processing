#!/bin/bash
#SBATCH -J APA_config
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

# Prepare config file
/data/taoyuhuan/projects/exOmics_RNA/bin/scripts/prepareDaParConfig.py -p ${output}/${positive}.txt -n ${output}/${negative}.txt -i ${output}/wig -o ${output}/dapar --config ${output}/config/config_${positive}vs${negative}.txt #python3 

# Run DaPar
#/data/taoyuhuan/tools/DaPars/dapars/src/DaPars_main.py ${output}/config_${positive}vs${negative}.txt

# Summarize result
#scripts/parseDapar.py -i ${output}/dapar/result_All_Prediction_Results.txt -c ${output}/config.txt -l ${output}/matrix/long.txt -s ${output}/matrix/short.txt -p ${output}/matrix/PDUI.txt
