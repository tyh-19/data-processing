#!/bin/bash
#SBATCH -J add_ID
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

input="/data/taoyuhuan/projects/PBMC_filtered_SNP/multiomics_paired"
output="/data/taoyuhuan/projects/PBMC_filtered_SNP/multiomics_paired_add_ID"
mkdir -p ${output}

echo "Start add ID for at `date`."
for sample in $(ls ${input} | grep rmPBMC | grep -v "CRC-PKU-5\|NC-PKU-mix16\|NC-PKU-mix17\|NC-PKU-mix19\|NC-PKU-mix20\|NC-PKU-mix22\|NC-PKU-mix30\|STAD-PKU-4" | cut -d "." -f 1)
do
echo "Start add ID ${sample} at `date`."
#gzip -d ${input}/${sample}.rmPBMC.SNP.gz
grep "^#" ${input}/${sample}.rmPBMC.SNP > ${output}/${sample}.rmPBMC.addID.SNP
grep -v "^#" ${input}/${sample}.rmPBMC.SNP | awk '{$3=$1"-"$2}1' OFS='\t' >> ${output}/${sample}.rmPBMC.addID.SNP
echo "End add ID ${sample} at `date`."
done

echo "Finish add ID at `date`."
