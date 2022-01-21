#!/bin/bash
#SBATCH -J Insert_length
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# path to /dataset/output/bam
input="/data/taoyuhuan/projects/exOmics_RNA/$1/output/bam"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Insert_length/$1"
mkdir ${output}
reference="/data/taoyuhuan/projects/exOmics_RNA/level_3_Insert_length/forMISO/exon_100"

#echo $(date +%F%n%T)
#exon_utils --get-const-exons /data/taoyuhuan/projects/exOmics_RNA/level_3_Insert_length/gencode.v34.annotation.gff3 --min-exon-size 100 --output-dir ${reference}

for file in $(ls ${input} | grep -v log | grep -v nohup)
    do
    echo ${file}
    mkdir ${output}/${file}
    echo $(date +%F%n%T)
    pe_utils --compute-insert-len ${input}/${file}/genome_rmdup.bam ${reference}/gencode.v34.annotation.min_100.const_exons.gff  --output-dir ${output}/${file}
done
