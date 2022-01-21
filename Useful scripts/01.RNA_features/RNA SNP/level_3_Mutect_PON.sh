#!/bin/bash
#SBATCH -J Mutect
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
#input="/data/taoyuhuan/projects/exOmics_RNA/${dataset}/output/bam"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}"
mkdir ${output}
mkdir ${output}/Mutect
mkdir ${output}/Mutect/PON
ref="/data/taoyuhuan/projects/exOmics_RNA/bin"
log="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/log"
mkdir ${log}
mkdir ${log}/Mutect

echo "Start build panel of normal for ${dataset} at `date`."
for sample in $2
do
echo "Call somatic mutation in ${sample}."
gatk Mutect2 -R ${ref}/genome/fasta/hg38.fa -I ${output}/SNP/bam_forSNP/${sample}_forSNP.bam -max-mnp-distance 0 -O ${output}/Mutect/PON/${sample}.vcf.gz
done
echo "Done at `date`."

#gatk GenomicsDBImport -R ${ref}/genome/fasta/hg38.fa --genomicsdb-workspace-path ${output}/Mutect/PON -V sample1.vcf.gz -V sample2.vcf.gz ...


