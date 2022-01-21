#local cibersort
source('Cibersort.R')
gene_expression_matrix <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/multiomics_paired_20211113_TPM.txt",sep = "\t")
gene_expression_matrix$`Gene.Length` <- as.character(lapply(strsplit(as.character(gene_expression_matrix$`Gene.Length`),".",fixed = TRUE),function(x) x[1]))

gene_expression_matrix <- aggregate(gene_expression_matrix[,-1], list(Gene=gene_expression_matrix[,1]), FUN = sum)
rownames(gene_expression_matrix) <- gene_expression_matrix$Gene
write.table(gene_expression_matrix[,-1],"TPM/multiomics_paired_20211113_TPM.txt",sep = "\t", quote = FALSE, row.names = TRUE,col.names = TRUE)

res_cibersort <- CIBERSORT("signature/Blood_CCLE_signature_genes_1_0.5_20.txt","TPM/multiomics_paired_20211113_TPM.txt",perm = 1000, QN = F)
write.table(res_cibersort, file="output/res_cibersort_multiomics_Blood_CCLE_1_0.5_20.txt", sep="\t", col.names=T, row.names=T, quote=F)
