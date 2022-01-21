#!/bin/bash
#SBATCH -J APA_wig
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
#positive=$2
#negative=$3
input="/data/taoyuhuan/projects/exOmics_RNA/${dataset}/output/bam"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_APA/${dataset}"
mkdir -p ${output}/wig
mkdir -p ${output}/dapar
mkdir -p ${output}/config
mkdir -p ${output}/matrix
mkdir -p ${output}/bam_sorted_by_coordinate

for sample in $(ls ${input})
do
echo "Sorting ${sample} at `date`:"
samtools sort -o ${output}/bam_sorted_by_coordinate/${sample}.bam -@ 16 ${input}/${sample}/genome_rmdup.bam
echo "Transcform to wig at `date`"
bedtools genomecov -ibam ${output}/bam_sorted_by_coordinate/${sample}.bam -bg -split | sort-bed - > ${output}/wig/${sample}.wig
done

# Prepare config file
#/data/taoyuhuan/projects/exOmics_RNA/bin/scripts/prepareDaParConfig.py -p ${output}/${positive}.txt -n ${output}/${negative}.txt -i ${output}/wig -o ${output}/dapar --config ${output}/config.txt #python3 

# Run DaPar
#/data/taoyuhuan/tools/DaPars/dapars/src/DaPars_main.py ${output}/config.txt

# Summarize result
#scripts/parseDapar.py -i ${output}/dapar/result_All_Prediction_Results.txt -c ${output}/config.txt -l ${output}/matrix/long.txt -s ${output}/matrix/short.txt -p ${output}/matrix/PDUI.txt
