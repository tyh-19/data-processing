#!/bin/bash
#SBATCH -J annotate_SNP
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

##Downlaod ensGene annotation (basic v31 from annovar)
##perl annotate_variation.pl -webfrom annovar -build hg38 -downdb ensGene humandb/hg38/

dataset=$1
software="/data/taoyuhuan/tools/annovar/annovar"
ref="/data/taoyuhuan/tools/annovar/annovar/humandb/hg38/"
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/SNP"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/SNP/annotation"
mkdir -p ${output}

echo "Start ${dataset} at `date`"

for sample in $(ls ${input} | grep "rmEDIT.SNP" |cut -d "." -f 1)
do
echo "Start ${sample} at `date`"
#vcf2avinput
perl ${software}/convert2annovar.pl -format vcf4 -withfreq ${input}/${sample}.rmEDIT.SNP > ${output}/${sample}.avinput

#annotation based on gene
perl ${software}/annotate_variation.pl -geneanno -dbtype ensGene -out ${output}/${sample} -build hg38 ${output}/${sample}.avinput ${ref}
done

echo "End at `date`"
