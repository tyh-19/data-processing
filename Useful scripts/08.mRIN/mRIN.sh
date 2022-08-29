#!/bin/bash
#SBATCH -J RNA_degradation
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

input_dir="/data/taoyuhuan/projects/exOmics_RNA/level_4/multiomics_paired/bedgraph"
work_dir="/data/taoyuhuan/projects/exOmics_RNA/level_4/test_mRIN/multiomics_all_genes"
RNA="long_RNA.gencode"

mkdir -p ${work_dir}/cdf
mkdir -p ${work_dir}/ks
mkdir -p ${work_dir}/ks_mat
mkdir -p ${work_dir}/res_mRIN

echo "#config file" > ${work_dir}/${RNA}_config.txt
for sample in $(ls ${input_dir} | grep -v "\.-\|+\|normalized" | grep "pico")
do
echo "${RNA}_${sample%*.bedgraph}.ks.txt	${sample%*.bedgraph}" >> ${work_dir}/${RNA}_config.txt
perl /data/taoyuhuan/tools/mRIN/gen_transcript_cdf.pl -v ${work_dir}/${RNA}.bed ${input_dir}/${sample} ${work_dir}/cdf/${RNA}_${sample%*.bedgraph}.cdf.bedgraph
perl /data/taoyuhuan/tools/mRIN/ks_test_uniform.pl -v ${work_dir}/cdf/${RNA}_${sample%*.bedgraph}.cdf.bedgraph ${work_dir}/ks/${RNA}_${sample%*.bedgraph}.ks.txt
rm ${work_dir}/cdf/${RNA}_${sample%*.bedgraph}.cdf.bedgraph
done

perl /data/taoyuhuan/tools/mRIN/gen_ks_matrix.pl -v -base ${work_dir}/ks --min-avg-cov 2 -v ${work_dir}/${RNA}_config.txt ${work_dir}/ks_mat/${RNA}_all_samples.KS.mat.txt

Rscript /data/taoyuhuan/tools/mRIN/cal_mrin.R -k ${work_dir}/ks_mat/${RNA}_all_samples.KS.mat.txt -m ${work_dir}/res_mRIN/${RNA}_out.mRIN.txt -G ${work_dir}/res_mRIN/${RNA}_out.GIS.txt -v
