#! /bin/bash
dataset=$1
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_TE/${dataset}/count_all"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_TE/${dataset}/count_all/count_matrix"
mkdir -p ${output}

first=$(ls ${input}|grep .cntTable| grep -v summary | head -n 1)
file=$(ls ${input}|grep cntTable| grep -v summary) 

echo ${first}
echo ${file}

#get gene id, delete first line(featurecount step generated)
cut -f 1 ${input}/${first} > ${output}/rowname.tmp.txt
cat ${output}/rowname.tmp.txt | tr -s '[:blank:]' '|' > ${output}/rowname.txt
#sed -i '1d' ${output}/rowname.txt
rm ${output}/rowname.tmp.txt
#get featurecount of each group, delete first line and alter second line to group info
for filename in ${file}
	do
	echo ${filename}
	cut -f 2 -d "	" ${input}/${filename} > ${output}/${filename}.tmp.txt
	sed -i '1d' ${output}/${filename}.tmp.txt
	#get the group info from filename. Since  cut only process file or dir, i touch a temp.txt to process filename.
	touch ${filename}.txt
	group=$(ls ${filename}.txt|cut -d "." -f 1)
	echo ${group}
	sed -i "1i \\${group}" ${output}/${filename}.tmp.txt
	rm ${filename}.txt
done

#merge to a matrix
cd ${output}
paste *.tmp.txt > Matrix_rawcount.txt

#merge rowname and matrix
paste rowname.txt  Matrix_rawcount.txt > TEtranscript_count_matrix

#convert to .csv
#echo "converting .txt to .csv"
#cat featurecount_summary.txt | tr -s '[:blank:]' ',' > featurecount_summary.csv 
#cat rowname.txt | tr -s '[:blank:]' ',' > rowname.csv
#cat Matrix_rawcount.txt | tr -s '[:blank:]' ',' > Matrix_rawcount.csv

#delete processing data
rm ${output}/*.txt
mv TEtranscript_count_matrix TEtranscript_count_matrix.txt

cat TEtranscript_count_matrix.txt | grep -v "ENSG00000198888\|ENSG00000198763\|ENSG00000198840\|ENSG00000212907\|ENSG00000198886\|ENSG00000198786\|ENSG00000198695\|ENSG00000198727\|ENSG00000198804\|ENSG00000198712\|ENSG00000198938\|ENSG00000198899\|ENSG00000228253\|ENSG00000211459\|ENSG00000210082\|ENSG00000210127\|ENSG00000210174\|ENSG00000210135\|ENSG00000210154\|ENSG00000210140\|ENSG00000210194\|ENSG00000210107\|ENSG00000210164\|ENSG00000210176\|ENSG00000210100\|ENSG00000209082\|ENSG00000210191\|ENSG00000210156\|ENSG00000210112\|ENSG00000210049\|ENSG00000210196\|ENSG00000210151\|ENSG00000210184\|ENSG00000210195\|ENSG00000210117\|ENSG00000210144\|ENSG00000210077" > TEtranscript_count_matrix_noMT.txt



