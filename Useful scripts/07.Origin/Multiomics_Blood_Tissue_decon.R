source('Cibersort.R')
res_cibersort <- CIBERSORT("signature/signature_genes_5_0.9_noCancer.txt","TPM/multiomics_paired_20211113_TPM.txt",perm = 1000, QN = F)
write.table(res_cibersort, file="output/res_cibersort_multiomics_noCancer.txt", sep="\t", col.names=T, row.names=T, quote=F)