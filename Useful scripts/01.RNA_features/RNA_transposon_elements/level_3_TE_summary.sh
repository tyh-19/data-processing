# !/bin/bash
dataset=$1
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_TE/${dataset}"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_TE/${dataset}/countmatrix"
mkdir -p ${output}

first=$(ls ${input} | grep -v "count\|log\|txt\|err\|out\|sh\|summary" | head -1)
cut -f 1 ${input}/${first}/${first}/quant.sf > ${output}/rownames.txt
sed -i "1iTE" ${output}/rownames.txt

for sample in $(ls ${input} | grep -v "count\|log\|txt\|err\|out\|sh\|summary")
do
echo ${sample}
cat ${input}/${sample}/${sample}/quant.sf | cut -f 4 > ${output}/${sample}.TPM.tmp.txt
sed -i "1i${sample}" ${output}/${sample}.TPM.tmp.txt

cat ${input}/${sample}/${sample}/quant.sf | cut -f 5 > ${output}/${sample}.readnum.tmp.txt
sed -i "1i${sample}" ${output}/${sample}.readnum.tmp.txt
done

paste ${output}/rownames.txt ${output}/*.TPM.tmp.txt > ${output}/${dataset}_TE_TPM.txt
paste ${output}/rownames.txt ${output}/*.readnum.tmp.txt > ${output}/${dataset}_TE_readnum.txt

rm ${output}/*.tmp.txt
rm ${output}/rownames.txt
