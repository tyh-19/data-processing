#!/bin/bash
#SBATCH -J ML_DIFF
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

workdir="/data/taoyuhuan/projects/Machine_learning/"
feature_matrix=$1
class="/data/taoyuhuan/projects/Machine_learning/20210722_multiomics/sample/class_pico.txt"
sample="/data/taoyuhuan/projects/Machine_learning/20210722_multiomics/sample/sample_ids_pico.txt"
positive=$2
negative=$3
outdir="/data/taoyuhuan/projects/Machine_learning/$4"

mkdir -p ${outdir}
mkdir -p ${outdir}/${positive}vs${negative}/

Rscript /data/taoyuhuan/projects/Machine_learning/scripts/differential_expression.R -i ${workdir}/${feature_matrix} -c ${class} -s ${sample} -p ${positive} -n ${negative} -m wilcox --norm-method NA --pseudo-count 1 -o ${outdir}/${positive}vs${negative}/Diff_wilcox.txt
