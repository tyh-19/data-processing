#!/bin/bash
#SBATCH -J DNA_SNP
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
tmp="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/tmp"
mkdir -p ${tmp}
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/SNP/bam_forSNP"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}"
mkdir -p ${output}/SNP
mkdir -p ${output}/SNP/tmp
ref="/data/taoyuhuan/projects/exOmics_RNA/bin"
log="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/log"
mkdir -p ${log}
mkdir -p ${log}/SNP

echo "Start DNA SNP analysis for ${dataset} at `date`."
for sample in $(ls ${input} | grep -v bai |grep -w $2)
do
echo "Start DNA SNP analysis ${sample} at `date`."
/data/taoyuhuan/tools/gatk-4.1.9.0/gatk HaplotypeCaller --java-options -Xmx4G \
        -R ${ref}/genome/fasta/hg38.fa \
        -I ${output}/SNP/bam_forSNP/${sample%.bam*}.bam -O ${output}/SNP/${sample%.bam*} \
        --tmp-dir ${output}/SNP/tmp 2>  ${log}/SNP/${sample%.bam*}.callSNP.log

/data/taoyuhuan/tools/gatk-4.1.9.0/gatk VariantFiltration --java-options -Xmx4G \
        -R ${ref}/genome/fasta/hg38.fa -V ${output}/SNP/${sample%.bam*} \
        -window 35 -cluster 3 \
        --tmp-dir ${tmp} \
        --filter-name FS20 -filter "FS > 20.0" \
        --filter-name QD2 -filter "QD < 2.0" \
        --filter-name DP10 -filter "DP < 10.0" \
        --filter-name QUAL20 -filter "QUAL < 20.0" -O ${output}/SNP/${sample%.bam*}.SNP 2> ${log}/SNP/${sample%.bam*}.filter.log
done
