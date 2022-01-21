#!/bin/bash
#SBATCH -J chimeric
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
strandness="no" # forward, reverse or no
input="/data/taoyuhuan/projects/exOmics_RNA/${dataset}/output/unmapped"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_chimeric/${dataset}"
mkdir ${output}
mkdir ${output}/counts

echo "Start counting ${dataset} chimeric RNA at `date`."
for sample in $(ls ${input})
do
echo "Mapping ${sample} chimeric RNA at `date`."
mkdir -p ${output}/${sample}
STAR --genomeDir /data/taoyuhuan/projects/exOmics_RNA/bin/genome/chimera-star-index \
          --readFilesIn ${input}/${sample}/circRNA_1.fastq.gz ${input}/${sample}/circRNA_2.fastq.gz \
          --runThreadN 16 \
          --outFileNamePrefix ${output}/${sample}/${sample} \
          --outSAMtype BAM Unsorted \
          --outReadsUnmapped Fastx \
          --readFilesCommand gzip -d -c \
          --outSAMmultNmax 1 \
          --seedPerWindowNmax 20

echo "Repairing ${sample} at `date`"
# To fix the bug of STAR_2.5.3a_modified 
/data/taoyuhuan/tools/bbmap/bbmap/repair.sh in=${output}/${sample}/${sample}Unmapped.out.mate1 in2=${output}/${sample}/${sample}Unmapped.out.mate2 out=${output}/${sample}/unmapped_1.fastq.gz out2=${output}/${sample}/unmapped_2.fastq.gz overwrite=t

echo "Counting ${sample} at `date`"
mkdir -p ${output}/counts/${sample}
/data/taoyuhuan/projects/exOmics_RNA/bin/scripts/count_reads.py count_circrna -s ${strandness} --paired-end -i ${output}/${sample}/${sample}Aligned.out.bam -o ${output}/counts/${sample}/${sample}


done
echo "Finished counting ${dataset} chimeric RNA at `date`."

