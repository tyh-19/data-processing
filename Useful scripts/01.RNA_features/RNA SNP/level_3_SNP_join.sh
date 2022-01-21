#!/bin/bash
#SBATCH -J annotate_SNP
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
software="/data/taoyuhuan/tools/annovar/annovar"
ref="/data/taoyuhuan/tools/annovar/annovar/humandb/hg38/"
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/SNP"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/SNP/matrix"
mkdir -p ${output}

cat | awk '{print $4,$10}' > tmp
sort 
join -a1 -a2 -o 1 2.1 test5.txt test6.txt | wc -l 
~
