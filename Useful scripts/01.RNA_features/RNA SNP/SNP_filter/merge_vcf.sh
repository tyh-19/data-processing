#!/bin/bash
#SBATCH -J merge_vcf
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

input="/data/taoyuhuan/projects/PBMC_filtered_SNP/PBMC_SNP"
output="/data/taoyuhuan/projects/PBMC_filtered_SNP"

first=$(ls ${input} | grep rmEDIT.sorted.vcf.gz | grep -v tbi | cut -d "." -f 1 | head -1)
second=$(ls ${input} | grep rmEDIT.sorted.vcf.gz | grep -v .tbi |  cut -d "." -f 1 | head -2 | tail -1)

echo "Merge ${first} and ${second} at `date`"
vcf-merge ${input}/${first}.rmEDIT.sorted.vcf.gz ${input}/${second}.rmEDIT.sorted.vcf.gz > ${output}/PBMC.vcf
vcf-sort ${output}/PBMC.vcf > ${output}/PBMC.sorted.vcf
bgzip -c ${output}/PBMC.sorted.vcf > ${output}/PBMC.sorted.vcf.gz
tabix -p vcf ${output}/PBMC.sorted.vcf.gz

for sample in $(ls ${input} | grep rmEDIT.sorted.vcf.gz | grep -v .tbi | cut -d "." -f 1 | grep -v ${first} | grep -v ${second})
do echo "${sample} at `date`"
vcf-merge ${input}/${sample}.rmEDIT.sorted.vcf.gz ${output}/PBMC.sorted.vcf.gz > ${output}/PBMC.tmp.vcf
rm ${output}/PBMC.vcf
mv ${output}/PBMC.tmp.vcf ${output}/PBMC.vcf

vcf-sort ${output}/PBMC.vcf > ${output}/PBMC.sorted.vcf
bgzip -c ${output}/PBMC.sorted.vcf > ${output}/PBMC.sorted.vcf.gz
tabix -p vcf ${output}/PBMC.sorted.vcf.gz

done
echo "Finish at `date`"
