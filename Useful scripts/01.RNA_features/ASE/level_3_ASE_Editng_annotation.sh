#!/bin/bash
#SBATCH -J annotate_ASE_Editing
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=%j.out
#SBATCH --error=%j.err

##Downlaod ensGene annotation (basic v31 from annovar)
##perl annotate_variation.pl -webfrom annovar -build hg38 -downdb ensGene humandb/hg38/

dataset=$1
target=$2
software="/data/taoyuhuan/tools/annovar/annovar"
ref="/data/taoyuhuan/tools/annovar/annovar/humandb/hg38/"

if [ ${target} = "COSMIC" ];then
echo "COSMIC start at `date`"

input="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/ASE/COSMIC"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/ASE/COSMIC/annotation"
mkdir -p ${output}

for sample in $(ls ${input} | grep "COSMIC" | cut -d "." -f 1)
do
echo "COSMIC Start ${sample} at `date`"
#gunzip
#mv ${input}/${sample}.COSMIC ${input}/${sample}.COSMIC.gz

#gzip -d ${input}/${sample}.COSMIC.gz

cat ${input}/${sample}.COSMIC | awk '{print $1,$2,$2,$4,$5}' | tr " " "	" > ${output}/${sample}.avinput

#annotation based on gene
perl ${software}/annotate_variation.pl -geneanno -dbtype ensGene -out ${output}/${sample} -build hg38 ${output}/${sample}.avinput ${ref}

done

echo "COSMIC End at `date`"

######################################

elif [ ${target} = "dbSNP" ]; then
echo "dbSNP Start at `date`"
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/ASE/dbSNP"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/ASE/dbSNP/annotation"
mkdir -p ${output}

echo "Start ${dataset} at `date`"

for sample in $(ls ${input} | grep "dbSNP" | cut -d "." -f 1 | head -247 | tail -2)
do
echo "dbSNP Start ${sample} at `date`"
#gunzip
#mv ${input}/${sample}.dbSNP ${input}/${sample}.dbSNP.gz

#gzip -d ${input}/${sample}.dbSNP.gz

cat ${input}/${sample}.dbSNP | awk '{print $1,$2,$2,$4,$5}' | tr " " "	" > ${output}/${sample}.avinput

#annotation based on gene
perl ${software}/annotate_variation.pl -geneanno -dbtype ensGene -out ${output}/${sample} -build hg38 ${output}/${sample}.avinput ${ref}

done

echo "dbSNP End at `date`"

######################################

elif [ ${target} = "REDI" ]; then
echo "REDIportal Start at `date`"
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/REDIportal"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/REDIportal/annotation"
mkdir -p ${output}

for sample in $(ls ${input} | grep "REDI" |cut -d "." -f 1)
do
echo "REDIportal Start ${sample} at `date`"
#gunzip
mv ${input}/${sample}.REDI ${input}/${sample}.REDI.gz

gzip -d ${input}/${sample}.REDI.gz

cat ${input}/${sample}.REDI | awk '{print $1,$2,$2,$4,$5}' | tr " " "	" > ${output}/${sample}.avinput

#annotation based on gene
perl ${software}/annotate_variation.pl -geneanno -dbtype ensGene -out ${output}/${sample} -build hg38 ${output}/${sample}.avinput ${ref}

done

echo "REDIportal End at `date`"

echo "End at `date`"

else
echo "Please pass 'COSMIC'/'dbSNP'/'REDI'." 
fi








