#!/bin/bash

input=$1
output=$1

echo "sample	human_percent	fungi_percent	bacteria_percent	archea_percent	virus_percent	human_readnum	fungi_readnum	bacteria_readnum	archea_readnum	virus_reaednum	genus_num specie_num" > ${output}/stat_kraken_report.txt

for sample in $(ls ${input}/report | grep -v nohup.out | grep -v stat_kraken_report.txt | grep -v stat_kraken_report_summary.txt)
do
echo ${sample}
human_percent=`cat ${input}/report/${sample} | awk -F " " '($5=="9606") {print $1}'` 
fungi_percent=`cat ${input}/report/${sample} | awk -F " " '($5=="4751") {print $1}'`
bacteria_percent=`cat ${input}/report/${sample} | awk -F " " '($5=="2") {print $1}'`
archea_percent=`cat ${input}/report/${sample} | awk -F " " '($5=="2157") {print $1}'`
virus_percent=`cat ${input}/report/${sample} | awk -F " " '($5=="10239") {print $1}'`

human_readnum=`cat ${input}/report/${sample} | awk -F " " '($5=="9606") {print $2}'`
fungi_readnum=`cat ${input}/report/${sample} | awk -F " " '($5=="4751") {print $2}'`
bacteria_readnum=`cat ${input}/report/${sample} | awk -F " " '($5=="2") {print $2}'`
archea_readnum=`cat ${input}/report/${sample} | awk -F " " '($5=="2157") {print $2}'`
virus_reaednum=`cat ${input}/report/${sample} | awk -F " " '($5=="10239") {print $2}'`

genus=`cat ${input}/report/${sample} | awk '($4=="G") {print $0}' | wc -l`
specie=`cat ${input}/report/${sample} | awk '($4=="S") {print $0}' | wc -l`

if [ -z ${human_percent} ]
	then
	human_percent=0
fi

if [ -z ${fungi_percent} ]
        then
        fungi_percent=0
fi

if [ -z ${bacteria_percent} ]
        then
        bacteria_percent=0
fi

if [ -z ${archea_percent} ]
        then
        archea_percent=0
fi

if [ -z ${virus_percent} ]
        then
        virus_percent=0
fi

if [ -z ${human_readnum} ]
        then
        human_readnum=0
fi

if [ -z ${fungi_readnum} ]
        then
        fungi_readnum=0
fi

if [ -z ${bacteria_readnum} ]
        then
        bacteria_readnum=0
fi

if [ -z ${archea_readnum} ]
        then
        archea_readnum=0
fi

if [ -z ${virus_readnum} ]
        then
        virus_readnum=0
fi

if [ -z ${genus} ]
        then
        genus=0
fi

if [ -z ${specie} ]
        then
        specie=0
fi

echo "${sample}	${human_percent}	${fungi_percent}	${bacteria_percent}	${archea_percent}	${virus_percent}	${human_readnum}	${fungi_readnum}	${bacteria_readnum}	${archea_readnum}	${virus_reaednum}	${genus}	${specie}" >> ${output}/stat_kraken_report.txt
done

cat ${output}/stat_kraken_report.txt | awk '{print $1,$2,$7,$3+$4+$5+$6,$8+$9+$10+$11,$12,$13 }' > ${output}/stat_kraken_report_summary.tmp

sed -i '1d' ${output}/stat_kraken_report_summary.tmp
sed -i '1i\sample	human_percent	human_read	microbe_percent	microbe_read	genus_num	specie_num' ${output}/stat_kraken_report_summary.tmp

tr " " "	" < ${output}/stat_kraken_report_summary.tmp > ${output}/stat_kraken_report_summary.txt
rm ${output}/stat_kraken_report_summary.tmp
