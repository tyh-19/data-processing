#!/bin/bash
for line in $(cat /data/taoyuhuan/projects/PBMC_filtered_SNP/multiomics_paired_GWAS/snp.txt);do echo ${line};sbatch search_SNP.sh ${line};sleep 5s;done
