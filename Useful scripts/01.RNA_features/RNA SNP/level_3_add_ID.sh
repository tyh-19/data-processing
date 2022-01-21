#!/bin/bash
#SBATCH -J filter_SNP
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/SNP"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}"
ref="/data/taoyuhuan/projects/exOmics_RNA/bin"

echo "Start add ID for ${dataset} at `date`."
for sample in $(ls ${input} | grep rmEDIT | grep -v "CRC-PKU-5\|NC-PKU-mix16\|NC-PKU-mix17\|NC-PKU-mix19\|NC-PKU-mix20\|NC-PKU-mix22\|NC-PKU-mix30\|STAD-PKU-4" | cut -d "." -f 1)
do
echo "Start add ID ${sample} at `date`."
#gzip -d ${input}/${sample}.rmEDIT.SNP
grep "^#" ${input}/${sample}.rmEDIT.SNP > ${output}/SNP/${sample}.rmEDIT.addID.SNP
grep -v "^#" ${input}/${sample}.rmEDIT.SNP | awk '{$3=$1"-"$2}1' OFS='\t' >> ${output}/SNP/${sample}.rmEDIT.addID.SNP
echo "End add ID ${sample} at `date`."
done

echo "Finish add ID for ${dataset} at `date`."
