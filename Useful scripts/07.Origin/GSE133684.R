#local cibersort
source('Cibersort.R')
gene_expression_matrix <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/GSE133684_exp_TPM-all.txt",sep = "\t",row.names = 1)
rownames(gene_expression_matrix) <- as.character(lapply(strsplit(rownames(gene_expression_matrix),"|",fixed = TRUE),function(x) x[1]))
write.table(gene_expression_matrix,"TPM/GSE133684_exp_TPM-all.txt",sep = "\t", quote = FALSE, row.names = TRUE,col.names = TRUE)

res_cibersort <- CIBERSORT("signature/blood_signature_genes.txt","TPM/GSE133684_exp_TPM-all.txt",perm = 1000, QN = F)
write.table(res_cibersort, file="output/res_cibersort_GSE133684_exp_TPM-all.txt", sep="\t", col.names=T, row.names=T, quote=F)