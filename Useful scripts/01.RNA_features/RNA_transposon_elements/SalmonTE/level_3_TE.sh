#!/bin/bash
#SBATCH -J TE
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
input="/data/taoyuhuan/projects/exOmics_RNA/${dataset}/output/cutadapt"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_TE/${dataset}"
mkdir ${output}

echo "TE for ${dataset}"

for sample in $(ls ${input} | grep 1.fastq.gz | cut -d "_" -f 1)
do
echo "Start analysing TE in ${sample} at `date`"
/data/taoyuhuan/projects/exOmics_RNA/bin/tools/SalmonTE/SalmonTE.py quant --reference=hs --num_threads=16 --outpath=${output}/${sample}  ${input}/${sample}_1.fastq.gz ${input}/${sample}_2.fastq.gz
echo "End analysing TE in ${sample} at `date`"
done

echo "${dataset} done at `date`"
