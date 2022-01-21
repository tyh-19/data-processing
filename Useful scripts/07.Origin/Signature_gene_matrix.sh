#!/bin/bash
#SBATCH -J 20_0.95_Signature_gene
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=%j.out
#SBATCH --error=%j.err

Rscript Signature_gene_matrix_v2.R -m ./Blood_CCLE/Blood_CCLE.txt -s ./Blood_CCLE/Blood_CCLE_sample_info.csv -e 1 -ts 0.9 -n 20 -o /data/taoyuhuan/projects/Origin/Blood_CCLE_titration

#Rscript Signature_gene_matrix.R -m ./Blood_GTEx_TCGA/Blueprint_TPM.sorted.txt -s ./Blood_GTEx_TCGA/Blueprint_sample_info.csv -e 5 -ts 0.8 -o /data/taoyuhuan/projects/Origin
