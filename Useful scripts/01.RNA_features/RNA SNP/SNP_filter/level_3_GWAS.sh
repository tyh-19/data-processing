# !/bin/bash
mkdir -p /data/taoyuhuan/projects/PBMC_filtered_SNP/multiomics_paired_GWAS
plink --bfile /data/taoyuhuan/projects/PBMC_filtered_SNP/multiomics_paired_bfile/summary/multiomics_paired --logistic --allow-extra-chr --adjust --out /data/taoyuhuan/projects/PBMC_filtered_SNP/multiomics_paired_GWAS/multiomics_paired_GWAS
