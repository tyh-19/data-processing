#!/bin/bash
#SBATCH -J annotate_SNP
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

##Downlaod ensGene annotation (basic v31 from annovar)
##perl annotate_variation.pl -webfrom annovar -build hg38 -downdb ensGene humandb/hg38/

software="/data/taoyuhuan/tools/annovar/annovar"
ref="/data/taoyuhuan/tools/annovar/annovar/humandb/hg38/"
input="/data/taoyuhuan/projects/PBMC_filtered_SNP/multiomics_paired"
output="/data/taoyuhuan/projects/PBMC_filtered_SNP/table_annovar"
mkdir -p ${output}

echo "Start annotating at `date`"

for sample in $(ls ${input} | grep "rmPBMC.SNP" | grep -v .gz |grep -v "CRC-PKU-mix1-pico\|NC-PKU-mix1-pico\|NC-PKU-mix2-pico\|NC-PKU-mix3-pico\|NC-PKU-mix4-pico\|NC-PKU-mix5-pico\|NC-PKU-mix6-pico\|NC-PKU-mix7-pico\|NC-PKU-mix8-pico\|NC-PKU-mix9-pico" |cut -d "." -f 1 | head -119 | tail -11)
do
echo "Start ${sample} at `date`"
#vcf2avinput
perl ${software}/convert2annovar.pl -format vcf4 -withfreq ${input}/${sample}.rmPBMC.SNP > ${output}/${sample}.avinput

#annotation based on gene
perl ${software}/table_annovar.pl ${output}/${sample}.avinput ${ref} -buildver hg38 -out ${output}/${sample} -remove -protocol ensGene -operation g -nastring . -thread 16
done

echo "End at `date`"
