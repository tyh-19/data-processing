# !/bin/bash

#input_dir is directory contains multiple samples' count results (use same gtf count by featurecounts)
input_dir=$1
log_dir_rownames=$2
log_dir_counts=$3
ref_genes=$4
output=$5

mkdir -p ${log_dir_rownames}
mkdir -p ${log_dir_counts}
mkdir -p ${output}

first=$(ls ${input_dir}|grep -v "log\|count" | grep -v summary | head -n 1)
all_countfiles=$(ls ${input_dir}|grep -v "log\|count" | grep -v summary)

echo "Using ${first}'s first column as matrix rowname."
echo "all samples in count matrix: ${all_countfiles}"

#get row names f r count matrix
cut -f 1,6 ${input_dir}/${first}/featurecount > ${log_dir_rownames}/rowname.tmp.txt
cat ${log_dir_rownames}/rowname.tmp.txt | tr -s '[:blank:]' '|' > ${log_dir_rownames}/rowname.txt
sed -i "1d" ${log_dir_rownames}/rowname.txt
rm ${log_dir_rownames}/rowname.tmp.txt

#add gene_type,gene_name to rowname
awk 'BEGIN { FS=OFS="|" } NR==FNR { hold[$1]=$0;next } {print $0, hold[$1] }' ${ref_genes} ${log_dir_rownames}/rowname.txt | awk 'BEGIN { FS=OFS="|" } {print $1,$5,$4,$2}' > ${log_dir_rownames}/annotated_rowname.txt
sed -i '1d' ${log_dir_rownames}/annotated_rowname.txt
sed -i '1igene_id|gene_symbol|gene_type|transcript_length' ${log_dir_rownames}/annotated_rowname.txt
sed -i 's/^|//' ${log_dir_rownames}/annotated_rowname.txt


#awk -F '|' 'NR==1;FNR==NR{a[$0];next}$1 in a{print $0}' OFS="|" ${log_dir_rownames}/rowname.txt  ${ref_genes} > ${log_dir_rownames}/annotated_rowname.txt

#join -1 1 -2 1 -t "|" ${log_dir_rownames}/rowname.txt ${ref_genes} > ${log_dir_rownames}/annotated_rowname.txt

#get all sample counts
for single_sample in ${all_countfiles}
do
echo ${single_sample}
cut -f 7 ${input_dir}/${single_sample}/featurecount > ${log_dir_counts}/${single_sample}.tmp.txt
sed -i '1,2d' ${log_dir_counts}/${single_sample}.tmp.txt
sed -i "1i${single_sample}" ${log_dir_counts}/${single_sample}.tmp.txt
done

#summarize matrix
paste ${log_dir_rownames}/annotated_rowname.txt ${log_dir_counts}/*.tmp.txt > ${output}/count_matrix.txt
