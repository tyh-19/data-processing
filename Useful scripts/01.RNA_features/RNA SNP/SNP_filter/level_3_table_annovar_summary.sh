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
input="/data/taoyuhuan/projects/PBMC_filtered_SNP/table_annovar"
output="/data/taoyuhuan/projects/PBMC_filtered_SNP/table_annovar"
mkdir -p ${output}

echo "Start ${dataset} at `date`"

for sample in $(ls ${output} | grep "avinput" |cut -d "." -f 1)
do
echo "Start ${sample} at `date`"
#vcf2avinput
#perl ${software}/convert2annovar.pl -format vcf4 -withfreq ${input}/${sample}.rmEDIT.SNP > ${output}/${sample}.avinput

#annotation based on gene
#perl ${software}/table_annovar.pl ${output}/${sample}.avinput ${ref} -buildver hg38 -out ${output}/${sample} -remove -protocol ensGene -operation g -nastring . -csvout -thread 16

cut -f '1-10' ${output}/${sample}.hg38_multianno.txt |sed '1d'|sed "s/$/\t${sample}/">> ${output}/all_sample.txt
done
sed -i '1s/^/Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tTumor_Sample_Barcode\n/' ${output}/all_sample.txt
echo "End at `date`"
