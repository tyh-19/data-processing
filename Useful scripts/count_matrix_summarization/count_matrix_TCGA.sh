# !/bin/bash

input_dir=$1
output_dir=$2
log_dir_rownames=$3
log_dir_counts=$4

mkdir -p ${output_dir}
mkdir -p ${log_dir_rownames}
mkdir -p ${log_dir_counts}

first=$(ls ${input_dir} | grep -v summary | head -n 1)
all_countfiles=$(ls ${input_dir} | grep -v summary | cut -d "." -f 1)

echo "Using ${first}'s first column as matrix rowname."
echo "all samples in count matrix: ${all_countfiles}"

#get row names f r count matrix
zcat ${input_dir}/${first} | cut -f 1 > ${log_dir_rownames}/rowname.txt
sed -i "1igene_id" ${log_dir_rownames}/rowname.txt

#get all sample counts
for single_sample in ${all_countfiles}
do
echo ${single_sample}
zcat ${input_dir}/${single_sample}.htseq.counts.gz | cut -f 2 > ${log_dir_counts}/${single_sample}.tmp.txt
sed -i "1i${single_sample}" ${log_dir_counts}/${single_sample}.tmp.txt
done

touch ${output_dir}/all_count.tmp.txt
for count in $(ls ${log_dir_counts} | grep -v "rowname")
do
paste ${output_dir}/all_count.tmp.txt ${log_dir_counts}/${count} > ${output_dir}/temp
cp ${output_dir}/temp ${output_dir}/all_count.tmp.txt
done

cut -f 2- -d "	" ${output_dir}/all_count.tmp.txt > ${output_dir}/all_count.txt
rm ${output_dir}/temp

paste ${log_dir_rownames}/rowname.txt ${output_dir}/all_count.txt > ${output_dir}/count_matrix.txt
