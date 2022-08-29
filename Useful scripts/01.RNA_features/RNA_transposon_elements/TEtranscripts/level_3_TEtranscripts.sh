#!/bin/bash
#SBATCH -J TEtranscripts
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

input="/data/taoyuhuan/projects/exOmics_RNA/multiomics_paired/output/bam"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_TE/multiomics_plasma"
mkdir -p ${output}/count_all
genic_GTF="/data/taoyuhuan/tools/exSEEK/genome/hg38/gtf/long_RNA.gtf"
TE_GTF="/data/taoyuhuan/tools/TEtranscripts-master/reference/hg38_rmsk_TE_20200804.gtf"
mkdir ${output}

echo "TE count for each sample"
for sample in $(ls ${input} | grep pico | head -137 | tail -17)
do
echo ${sample}
TEcount -b ${input}/${sample}/genome.bam \
              --GTF ${genic_GTF} \
              --TE  ${TE_GTF} \
              --format BAM \
              --stranded reverse \
              --project multiomics_plasma/count_all/${sample} \
              --mode multi
echo "${sample} finished at `date`"
done

#echo "TEtranscripts test by DESeq2"
#TEtranscripts -t ${input}/STAD-PKU-8-pico/genome_rmdup.bam ${input}/STAD-PKU-9-pico/genome_rmdup.bam\
#              -c ${input}/NC-PKU-10-pico/genome_rmdup.bam ${input}/NC-PKU-11-pico/genome_rmdup.bam\
#              --GTF ${genic_GTF} \
#              --TE  ${TE_GTF} \
#              --format BAM \
#              --stranded reverse \
#              --project multiomics_plasma/test \
#              --mode multi
#echo "Jobs done at `date`"
