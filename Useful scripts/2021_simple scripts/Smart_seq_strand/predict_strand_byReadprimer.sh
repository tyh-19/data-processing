#!/bin/bash
#SBATCH -J strand
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --output=%j.out
#SBATCH --error=%j.err


# This script only works for SMARTer Stranded Total RNA-Seq Kit—Pico Input Mammalian (to distinguish V1 or V2)
Rd1_c="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
Rd2_c="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
input="/data/taoyuhuan/projects/pico_strand/test"
tmp="/data/taoyuhuan/projects/pico_strand/tmp"

for sample in $(ls ${input} | grep "1.fastq.gz" | cut -d "_" -f 1)
do
echo ${sample} "at `date`"
zcat ${input}/${sample}_1.fastq.gz | grep ${Rd2_c} | head -1000 | while read line
do
read1=${line%AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC*}
echo ${read1:0-3} >> ${tmp}/${sample}_1.tmp
done

sed -i '/^\s*$/d' ${tmp}/${sample}_1.tmp
file1_CCC=$(cat ${tmp}/${sample}_1.tmp | head -100 | grep CCC | wc -l)
echo ${file1_CCC} "in R1 first 1000 read through"

if ((${file1_CCC}>=13));then
	file1_predict="Reverse"
	else
	file1_predict="Forward"
fi

echo ${file1_predict}
echo "Finish 50%."

zcat ${input}/${sample}_2.fastq.gz | grep ${Rd1_c} | head -1000 | while read line
do
read2=${line%AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*}
echo ${read2:0-3} >> ${tmp}/${sample}_2.tmp
done

sed -i '/^\s*$/d' ${tmp}/${sample}_2.tmp
file2_CCC=$(cat ${tmp}/${sample}_2.tmp | head -100 | grep CCC | wc -l)
echo ${file2_CCC} "in R2 first 1000 read through"

if ((${file2_CCC}>=13));then
        file2_predict="Forward"
        else
        file2_predict="Reverse"
fi

echo ${file2_predict}

if (("${file1_predict}"=="${file2_predict}"));then
	echo -n ${sample} ":"
	predict=${file1_predict}
	else
	predict="Uncertain"
fi
echo "Finish 100%."
echo ${predict}
echo ${sample}":"${predict} >> ${tmp}/strand.txt 
done
