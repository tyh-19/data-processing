#!/bin/bash
#SBATCH -J filter_PBMC
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

input="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/multiomics_paired/SNP"
output="/data/taoyuhuan/projects/PBMC_filtered_SNP/multiomics_paired"
mkdir -p ${output}
ref="/data/taoyuhuan/projects/PBMC_filtered_SNP"

for sample in $(ls ${input} | grep rmEDIT | grep -v addID | grep "mix.-.-" | cut -d "." -f 1)
do
bedtools intersect -v -header -a ${input}/${sample}.rmEDIT.SNP -b ${ref}/PBMC.sorted.vcf.gz | bgzip -c > ${output}/${sample}.rmPBMC.SNP.gz
done
