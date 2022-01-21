#local cibersort
source('Cibersort.R')
gene_expression_matrix <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/GSE142987_sample_count_matrix_TPM.txt",sep = "\t")
rownames(gene_expression_matrix) <- as.character(lapply(strsplit(rownames(gene_expression_matrix),"|",fixed = TRUE),function(x) x[1]))
write.table(gene_expression_matrix,"TPM/GSE142987_sample_count_matrix_TPM.txt",sep = "\t", quote = FALSE, row.names = TRUE,col.names = TRUE)

res_cibersort <- CIBERSORT("signature/blood_signature_genes.txt","TPM/GSE142987_sample_count_matrix_TPM.txt",perm = 1000, QN = F)
write.table(res_cibersort, file="output/res_cibersort_GSE142987_sample_count_matrix_TPM.txt", sep="\t", col.names=T, row.names=T, quote=F)