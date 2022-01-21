# !/bin/bash

dataset=$1
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/ASE/COSMIC/annotation"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/ASE/COSMIC/matrix"
mkdir -p ${output}
mkdir -p ${output}/tmp
mkdir -p ${output}/all_SNP
#get all SNPs
cat ${input}/*.exonic_variant_function | grep -v unknown | grep "stopgain\|stoploss" | awk '{print $3,$2,$2,$4,$5,$6,$7,$8}' | sort | uniq > ${output}/all_SNP/all_SNP_stop.txt
cat ${input}/*.exonic_variant_function | grep -v unknown | grep -v "stopgain\|stoploss" | awk '{print $4,$2,$3,$5,$6,$7,$8,$9}' | sort | uniq > ${output}/all_SNP/all_SNP_nostop.txt
rm ${output}/all_SNP/all_SNP_${dataset}.txt #in case of sed -i change the file
cat ${output}/all_SNP/all_SNP_nostop.txt ${output}/all_SNP/all_SNP_stop.txt | sort | tr " " "	" > ${output}/all_SNP/all_SNP_${dataset}.txt
sed -i "1igene	consequence	variation	chromosome	start	end	reference	alternative" ${output}/all_SNP/all_SNP_${dataset}.txt

#different SNP types
#frameshift deletion
#frameshift insertion
#nonframeshift deletion
#nonframeshift insertion
#nonsynonymous SNV
#stopgain stopgain
#stoploss stoploss
#synonymous SNV

#get RP-RNA SNP
#cat ${output}/all_SNP/all_SNP_${dataset}.txt | grep "\<RPL\|\<RPS" | sort > ${output}/all_SNP/subset_SNP_${dataset}.txt

#process sample annotation
for sample in $(ls ${input} | grep "exonic_variant_function" | cut -d "." -f 1)
do
echo "Start ${sample} at `date`"
cat ${input}/${sample}.exonic_variant_function | grep -v unknown | grep -v stop | awk '{print $4,$10}' | tr " " "	" | sort > ${output}/tmp/${sample}_nostop.txt
cat ${input}/${sample}.exonic_variant_function | grep -v unknown | grep "stop" | awk '{print $3,$9}' | tr " " "	" | sort > ${output}/tmp/${sample}_stop.txt
cat ${output}/tmp/${sample}_stop.txt ${output}/tmp/${sample}_nostop.txt | sort > ${output}/tmp/${sample}.txt

#join
join -a1 -1 1 -2 1 -o 1.1 2.2 ${output}/all_SNP/subset_SNP_${dataset}.txt ${output}/tmp/${sample}.txt | tr " " "	" > ${output}/tmp/${sample}_in_all.txt
sed -i "1iSNP	${sample}" ${output}/tmp/${sample}_in_all.txt

cat ${output}/tmp/${sample}_in_all.txt | awk '{print $2}' > ${output}/tmp/${sample}_AF_in_all.txt
done

#simplified matrix, of which row index is gene symbols(with repeats). 
cat ${output}/all_SNP/subset_SNP_${dataset}.txt | awk -F ":" '{print $1}' > ${output}/tmp/1_rownames.txt
sed -i "1iSNV" ${output}/tmp/1_rownames.txt
paste ${output}/tmp/1_rownames.txt ${output}/tmp/*_AF_in_all.txt | tr " " "	" > ${output}/SNP_matrix.txt

#detailed matrix
sed -i "1igene  consequence     variation       chromosome      start   end     reference       alternative" ${output}/all_SNP/subset_SNP_${dataset}.txt
paste ${output}/all_SNP/subset_SNP_${dataset}.txt ${output}/tmp/*_AF_in_all.txt | tr " " "	" > ${output}/SNP_matrix_with_annotation.txt
echo "Done at `date`"

