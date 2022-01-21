#!/bin/bash
#SBATCH -J dapar
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

echo "Start calculate APA for ${dataset} between positive:${positive} and negative:${negative}"
echo "`date`"
# Run DaPar
python /data/taoyuhuan/tools/DaPars/dapars/src/DaPars_main.py ${output}/config/config_${positive}vs${negative}.txt

echo "End at `date`"
