# !/bin/bash

ref="/data/taoyuhuan/projects/exOmics_RNA/multiomics_paired/scripts/QC_scripts/bed/gene.bed"
bed="/data/taoyuhuan/projects/Coverage/bed"
gene="/data/taoyuhuan/projects/Coverage/genes/COSIMIC_CRC.txt"
bam="/data/taoyuhuan/projects/Coverage/CRC-PKU-32/bam"
output="/data/taoyuhuan/projects/Coverage/CRC-PKU-32"

cat ${gene} | head -4 | while read line
do
echo ${line}
mkdir -p ${output}/${line}
	for sample in $(ls ${bam} | grep ".bam" | cut -d "." -f 1)
	do
	echo ${sample}
	cat ${ref} | grep ${line} > ${bed}/${line}.bed
	bedtools intersect -a ${bam}/${sample}.bam -b ${bed}/${line}.bed | samtools sort > ${output}/${line}/${sample}_sorted.bam
	samtools depth ${output}/${line}/${sample}_sorted.bam -@ 12 -o ${output}/${line}/${sample}_sorted.txt

	total_mapped=$[$(samtools view -c -f 0x2 ${bam}/${sample}.bam)/2]
	
	echo ${total_mapped} > ${output}/${line}/${sample}_mapped_read_pairs.txt
	echo "${sample} mapped read-pairs:${total_mapped}"
	cat ${output}/${line}/${sample}_sorted.txt | awk -v value=${total_mapped} '{print $1,$2,$3/value*1000000}' | tr " " "\t" > ${output}/${line}/${sample}_sorted_cpm.txt
	done
done


