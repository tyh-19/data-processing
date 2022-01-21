# !/bin/bash

dataset=$1
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_chimeric/${dataset}/counts"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_chimeric/${dataset}/count_matrix"
mkdir -p ${output}

first=$(ls ${input} | grep -v "count\|log\|txt\|err\|out\|sh\|summary" | head -1)
cut -f 1 ${input}/${first}/${first} > ${output}/rownames.txt
sed -i "1iChimeric" ${output}/rownames.txt

for sample in $(ls ${input} | grep -v "count\|log\|txt\|err\|out\|sh\|summary")
do
echo ${sample}
cat ${input}/${sample}/${sample} | cut -f 2 > ${output}/${sample}.tmp.txt
sed -i "1i${sample}" ${output}/${sample}.tmp.txt
done

paste ${output}/rownames.txt ${output}/*.tmp.txt > ${output}/${dataset}_chimeric.txt

rm ${output}/*.tmp.txt
rm ${output}/rownames.txt
