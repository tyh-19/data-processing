#!/bin/bash
#SBATCH -J Combination
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

input_dir="/data/taoyuhuan/projects/Machine_learning/20211004_top50"
output_dir="/data/taoyuhuan/projects/Machine_learning/20211004_top50/Merged_seed6"
mkdir -p ${output_dir}
positive="STAD"
negative="NC"
Alteration_number=$1

Rscript /data/taoyuhuan/projects/Machine_learning/scripts/Combination.R -i ${input_dir} -f $1 -p ${positive} -n ${negative} -o ${output_dir} > ${output_dir}/test_MergeModel.log
