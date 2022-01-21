#!/bin/bash
#SBATCH -J Coverage
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=%j.out
#SBATCH --error=%j.err

ref="/data/taoyuhuan/projects/Coverage/gene.bed"
bed="/data/taoyuhuan/projects/Coverage/bed"
gene="/data/taoyuhuan/projects/Coverage/genes/COSIMIC_CRC.txt"
bam="/data/taoyuhuan/projects/Coverage/STAD-PKU-13/bam"
output="/data/taoyuhuan/projects/Coverage/STAD-PKU-13"

mkdir -p ${bed}

cat ${gene} | grep ENSG | while read line
do
echo ${line}
mkdir -p ${output}/${line}
	for sample in $(ls ${bam} | grep "qia.bam" | cut -d "." -f 1)
	do
	echo ${sample}
	cat ${ref} | grep ${line} > ${bed}/${line}.bed
	bedtools intersect -a ${bam}/${sample}.bam -b ${bed}/${line}.bed | samtools sort > ${output}/${line}/${sample}_sorted.bam
	samtools depth ${output}/${line}/${sample}_sorted.bam -o ${output}/${line}/${sample}_sorted.txt
	if [[ "$(echo ${sample} | grep "qia\|small")" != "" ]];then
		total_mapped=$[$(samtools view -c ${bam}/${sample}.bam)]
	else
		total_mapped=$[$(samtools view -c -f 0x2 ${bam}/${sample}.bam)/2]
	fi
	echo ${total_mapped} > ${output}/${line}/${sample}_mapped_read_pairs.txt
	echo "${sample} mapped read-pairs:${total_mapped}"
	cat ${output}/${line}/${sample}_sorted.txt | awk -v value=${total_mapped} '{print $1,$2,$3/value*1000000}' | tr " " "\t" > ${output}/${line}/${sample}_sorted_cpm.txt
	done
done


