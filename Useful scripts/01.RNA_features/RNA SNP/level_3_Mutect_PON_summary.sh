#!/bin/bash
#SBATCH -J Mutect
#SBATCH -p CN_BIOT
#SBATCH --nodes=5
#SBATCH --ntasks=80
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
#input="/data/taoyuhuan/projects/exOmics_RNA/${dataset}/output/bam"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}"
mkdir ${output}
mkdir ${output}/Mutect
mkdir ${output}/Mutect/PON
mkdir ${output}/Mutect/pon_db
ref="/data/taoyuhuan/projects/exOmics_RNA/bin"
log="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/log"
mkdir ${log}
mkdir ${log}/Mutect

echo "Start build panel of normal for ${dataset} at `date`."
gatk GenomicsDBImport -R ${ref}/genome/fasta/hg38.fa \
--genomicsdb-workspace-path pon_db \
-L /data/taoyuhuan/tools/exSEEK/genome/hg38/bed/genome.bed \
-V ${output}/Mutect/PON/CRC-2384058-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-2399129-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-27-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-28-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-29-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-30-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-32-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-34-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-35-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-36-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-37-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-38-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-39-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-40-N.vcf.gz \
-V ${output}/Mutect/PON/CRC-PKU-41-N.vcf.gz 

gatk CreateSomaticPanelOfNormals \
-R ${ref}/genome/fasta/hg38.fa \
-V gendb://pon_db \
-O ${output}/Mutect/pon.vcf.gz
echo "Done at `date`."

