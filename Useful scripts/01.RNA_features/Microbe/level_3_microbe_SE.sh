#!/bin/bash
#SBATCH -J SE_microbe
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err


# copy unmaped fastq.gz to ${input}/{dataset}/out/unmapped

dataset=$1
input="/data/taoyuhuan/projects/exOmics_RNA/${dataset}/output/unmapped"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_microbe/${dataset}"
mkdir ${output}
mkdir ${output}/unclassified
mkdir ${output}/report
database="/data/taoyuhuan/reference/kraken2-standard-db"

echo "Start counting ${dataset} microbe abundance at `date`."
for sample in $(ls ${input})
do
echo "Counting ${sample} microbe abundance at `date`."
kraken2 --db ${database} \
	--threads 16 \
	--unclassified-out ${output}/unclassified/${sample}_unclassified.fq \
	--report ${output}/report/${sample} \
	--use-names \
        ${input}/${sample}/circRNA.fa.gz
echo "Finished counting ${sample} microbe abundance at `date`."
done
echo "Finished counting ${dataset} microbe abundance at `date`."
