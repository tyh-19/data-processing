# !/bin/bash

input=/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/$1/SNP
output=/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/$1/SNP/summary
rm -r ${output}
mkdir ${output}

for sample in $(ls ${input} | grep addID | cut -d "." -f 1)
do
echo ${sample}

#--make-bed means to bed/bim/fam file
plink --vcf ${input}/${sample}.rmEDIT.addID.SNP --maf 0.05 --double-id --allow-extra-chr --make-bed --out ${output}/${sample}_tmp
# --recode means to ped/map file
#plink --vcf ${input}/${sample}.rmEDIT.SNP --maf 0.05 --double-id --allow-extra-chr --recode --out ${output}/${sample}

#exclude 3+ alleles
#plink --bfile ${output}/${sample}_tmp --exclude ${output}/${sample}.missnp --allow-extra-chr --make-bed --out ${output}/${sample}


echo "${output}/${sample}_tmp" >> ${output}/all_sample_tmp.txt

done

sed -i "1d" ${output}/all_sample_tmp.txt

# binary file
plink --noweb --bfile ${output}/$(ls ${input} | grep addID | cut -d "." -f 1|head -1)_tmp --merge-list ${output}/all_sample_tmp.txt --allow-extra-chr --out ${output}/$1_tmp
# ped/map
#plink --noweb --file ${output}/$(ls ${input} | grep rmEDIT | cut -d "." -f 1|head -1) --merge-list ${output}/all_sample.txt --allow-extra-chr --out ${output}/$1

##rm 3+ alleles
for sample in $(ls ${input} | grep addID | cut -d "." -f 1)
do
echo ${sample}

#exclude 3+ alleles
plink --bfile ${output}/${sample}_tmp --exclude ${output}/$1_tmp.missnp --allow-extra-chr --make-bed --out ${output}/${sample}
echo "${output}/${sample}" >> ${output}/all_sample.txt

done

sed -i "1d" ${output}/all_sample.txt

# binary file
plink --noweb --bfile ${output}/$(ls ${input} | grep addID | cut -d "." -f 1|head -1) --merge-list ${output}/all_sample.txt --allow-extra-chr --out ${output}/$1


