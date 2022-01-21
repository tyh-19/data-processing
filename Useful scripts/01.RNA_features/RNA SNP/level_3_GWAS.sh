# !/bin/bash
plink --bfile $1/SNP/summary/$1 --logistic --allow-extra-chr --adjust --out $1/SNP/summary/$1_GWAS
