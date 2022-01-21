#!/bin/bash

input="/data/taoyuhuan/projects/exOmics_RNA/level_3_Insert_length/$1"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Insert_length/$1/summary"
mkdir ${output}
tmp="/data/taoyuhuan/projects/exOmics_RNA/level_3_Insert_length/$1/tmp"
mkdir ${tmp}

file=$(ls ${input} | grep -v .err | grep -v .out | grep -v tmp | grep -v summary)

echo ${file}
for filename in ${file}
do
echo ${filename}
cat ${input}/${filename}/genome_rmdup.bam.insert_len | awk -F '	' '{print $2}' | awk 'BEGIN{RS=","}{print $0}' > ${tmp}/${filename}.tmp
sed -i "1,2d" ${tmp}/${filename}.tmp
sort -n ${tmp}/${filename}.tmp > ${tmp}/${filename}.sorted.tmp
sed -i "1i${filename}" ${tmp}/${filename}.sorted.tmp
done

paste ${tmp}/*.sorted.tmp > ${output}/miso_Insert_length_summary.txt
rm ${tmp}/*.tmp
