#!/bin/bash
#SBATCH -J search_SNP
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

input="/data/taoyuhuan/projects/PBMC_filtered_SNP/multiomics_paired_add_ID"
target=$1
output="/data/taoyuhuan/projects/PBMC_filtered_SNP/classification"
mkdir -p ${output}

echo ${target} > ${output}/${target}.txt
for sample in $(ls ${input} | grep rmPBMC | grep -v "CRC-PKU-5\|NC-PKU-mix16\|NC-PKU-mix17\|NC-PKU-mix19\|NC-PKU-mix20\|NC-PKU-mix22\|NC-PKU-mix30\|STAD-PKU-4" | cut -d "." -f 1)
do
echo -e ${sample}"\t\c" >> ${output}/${target}.txt
exist=$(less ${input}/${sample}.rmPBMC.addID.SNP | grep "${target}" | cut -d "	" -f 3)
AF=$(less ${input}/${sample}.rmPBMC.addID.SNP | grep "${target}" | cut -d " " -f 8 | cut -d ";" -f 2 | cut -d "=" -f 2)
if [ "${exist}" = "${target}" ];then
	echo -e ${AF}"\t\c" >> ${output}/${target}.txt
	echo ${AF} | awk '{if($0>=0.05&&$0<1) print "Exist!";else print "Not exist because MAF < 0.05 or AF=1"}' >> ${output}/${target}.txt
else
	echo -e "NA\t\c" >> ${output}/${target}.txt
	echo "Not Exist." >> ${output}/${target}.txt
fi
done
