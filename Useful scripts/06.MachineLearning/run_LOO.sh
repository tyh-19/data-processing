#!/bin/bash
#SBATCH -J ML_LOO
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

matrix="/data/taoyuhuan/projects/Machine_learning/$1"
featureNumber=$2
positive=$3
negative=$4
outdir=$5

mkdir -p ${outdir}
mkdir -p ${outdir}/${positive}vs${negative}/

Rscript /data/taoyuhuan/projects/Machine_learning/scripts/LOO_v3_yuhuan.R -i ${matrix} -f ${featureNumber} -p ${positive} -n ${negative} -o ${outdir}/${positive}vs${negative}/
