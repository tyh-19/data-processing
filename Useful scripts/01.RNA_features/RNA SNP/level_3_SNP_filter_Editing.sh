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

echo "Start filter Editing from SNP for ${dataset} at `date`."
for sample in $(ls ${input} | grep SNP | grep -v idx | grep -v bam_forSNP | cut -d "." -f 1)
do
echo "Start filter ${sample} at `date`."
## -v means only keep variants not appear in range specified in -b file
bedtools intersect -v -header -a ${output}/SNP/${sample}.SNP -b ${ref}/genome/vcf/REDIportal.vcf.gz | bgzip -c > ${output}/SNP/${sample}.rmEDIT.SNP.gz

echo "End filter ${sample} at `date`."
done

echo "Finish Editing ASE SNP analysis for ${dataset} at `date`."
