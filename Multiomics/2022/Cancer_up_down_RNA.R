library(edgeR)
#commen function for differential analysis: wilcox test
Wilcox_test <- function(mat_raw,des,output_res){
  norm_method <- 'NA'
  samples <- des$samples
  group <- des$group
  positive <- as.character(des[which(des$group=="positive"),]$samples)
  negative <- as.character(des[which(des$group=="negative"),]$samples)
  
  if(norm_method == 'NA' ){
    message('Matrix output without normalization.')
    mat_raw[is.na(mat_raw)] <- 0
    matrix <- mat_raw
  }else{
    message('Matrix output normalized by:',norm_method)
    matrix <- cpm(mat_raw, method=norm_method)
  }
  
  matrix <- matrix[,as.character(samples)]
  
  test_func <- function(x){
    positive_mean<-mean(x[group=="positive"])
    negative_mean<-mean(x[group=="negative"])
    if (positive_mean == negative_mean) {
      wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='two.sided')$p.value
    } else if (positive_mean > negative_mean) {
      wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='greater')$p.value
    } else if (positive_mean < negative_mean) {
      wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='less')$p.value
    }
  }
  
  #filter events < 20% samples in minial group
  #positive_prop <- rowSums(matrix[,positive] > 0)/length(positive)
  #negative_prop <- rowSums(matrix[,negative] > 0)/length(negative)
  #matrix <- matrix[(negative_prop>=0.2)&(positive_prop>=0.2),]
  
  #write.table(matrix, output_matrix, sep='\t', quote=FALSE, row.names=TRUE)
  
  pvalues <- apply(matrix, 1, test_func)
  matrix_logcpm = log2(matrix + 1)
  logFC <- apply(matrix_logcpm[,positive], 1, mean) -
    apply(matrix_logcpm[,negative], 1, mean)
  deltaFraction <- apply(matrix[,positive], 1, mean) -
    apply(matrix[,negative], 1, mean)
  
  positive_gini <- as.numeric(gini(t(matrix[,positive])))
  negative_gini <- as.numeric(gini(t(matrix[,negative])))
  
  res <- data.frame(log2FoldChange=logFC,
                    deltaFraction=deltaFraction,
                    pvalue=pvalues, 
                    padj=p.adjust(pvalues, method='BH'),
                    baseMean=apply(matrix, 1, mean),
                    positive_gini = positive_gini,
                    negative_gini = negative_gini)
  res <- res[order(res$pvalue,decreasing = FALSE),]
  write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE)
  res
}

setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx")

#CCLE
{
mat_raw <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/CCLE/CCLE_RNAseq_rsem_genes_tpm_20180929.txt",sep = "\t",header = TRUE,check.names = FALSE)
des <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/CCLE/Colon_cell_info.csv")
#des <- data.frame(apply(des,2,as.character)) 
output_res <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/CCLE/Colon_cell_DE.txt"
output <- Wilcox_test(mat_raw,des,output_res)
}

#TCGA+GTEx data integration
{
  #data matrix
  {
  mat_raw_COAD <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/TCGA/COAD/TCGA_COAD_count_matrix_TPM.txt",sep = "\t",header = TRUE,check.names = FALSE)
  mat_raw_READ <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/TCGA/READ/TCGA_READ_count_matrix_TPM.txt",sep = "\t",header = TRUE,check.names = FALSE)
  mat_raw_STAD <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/TCGA/STAD/TCGA_STAD_count_matrix_TPM.txt",sep = "\t",header = TRUE,check.names = FALSE)
  mat_raw_TCGA <- cbind(mat_raw_COAD,mat_raw_READ,mat_raw_STAD)
  mat_raw_TCGA$ensembl_gene_id <- as.character(lapply(strsplit(rownames(mat_raw_TCGA),"|",fixed = TRUE),function(x) x[1]))
  mat_raw_GTEx <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/GTEx/GTEx_ensembl_gene_tpm.txt",sep = "\t",header = TRUE,check.names = FALSE)
  mat_raw_GTEx$ensembl_gene_id <- lapply(strsplit(as.character(mat_raw_GTEx$Name),".",fixed = TRUE),function(x) x[1])
  }
  
  #sample info
  {
  des_COAD <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/TCGA/COAD/TCGA_COAD_PrimaryTumor_vs_TissueNormal.csv")
  des_READ <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/TCGA/READ/TCGA_READ_PrimaryTumor_vs_TissueNormal.csv")
  des_STAD <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/TCGA/STAD/TCGA_STAD_PrimaryTumor_vs_TissueNormal.csv")
  des_GTEx <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/GTEx/GTEx_noBlood_sample_info.csv")
  des_GTEx_Stomach <- des_GTEx[which(des_GTEx$group=="Stomach"),]
  des_GTEx_Colon <- des_GTEx[which(des_GTEx$group=="Colon"),]
  des_GTEx <- rbind(des_GTEx_Colon,des_GTEx_Stomach)
  }
  
  #remove unmatched samples
  unmatch_samples <- setdiff(as.character(des_GTEx$sample),colnames(mat_raw_GTEx))
  des_GTEx <- des_GTEx[!des_GTEx$sample %in% unmatch_samples,]
  write.table(as.data.frame(unmatch_samples) ,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/GTEx_unmatched_samples.txt",sep = "\t",quote=FALSE)
  
  #select samples
  mat_raw_GTEx_selected <- mat_raw_GTEx[,as.character(des_GTEx$sample),]
  rownames(mat_raw_GTEx_selected) <- mat_raw_GTEx$Name
  mat_raw_GTEx_selected$ensembl_gene_id <- as.character(lapply(strsplit(rownames(mat_raw_GTEx_selected),".",fixed = TRUE),function(x) x[1]))
  #merge TCGA and GTEx
  mat_raw_TCGA_GTEx <- left_join(mat_raw_TCGA,mat_raw_GTEx_selected,by=c("ensembl_gene_id"="ensembl_gene_id"))
  mat_raw_TCGA_GTEx[is.na(mat_raw_TCGA_GTEx)] <- 0
  
  mat_raw_TCGA_GTEx <- aggregate(mat_raw_TCGA_GTEx[,-which(colnames(mat_raw_TCGA_GTEx)=="ensembl_gene_id")], list(Gene=as.character(mat_raw_TCGA_GTEx$ensembl_gene_id)), FUN = sum)
  
  rownames(mat_raw_TCGA_GTEx) <- mat_raw_TCGA_GTEx$Gene
  mat_raw_TCGA_GTEx <- mat_raw_TCGA_GTEx[,-which(colnames(mat_raw_TCGA_GTEx)=="Gene")]
  write.table(mat_raw_TCGA_GTEx,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/mat_TPM_TCGA_GTEx.txt",sep = "\t",quote=FALSE)
}

#Colorectal cancer vs. normal tissue
{
  des_COAD <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/TCGA/COAD/TCGA_COAD_PrimaryTumor_vs_TissueNormal.csv")
  des_READ <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/TCGA/READ/TCGA_READ_PrimaryTumor_vs_TissueNormal.csv")
  des_TCGA_colorectal <- rbind(des_COAD,des_READ)
  des_TCGA_colorectal$batch <- "TCGA"
  des_GTEx <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/GTEx/GTEx_noBlood_sample_info.csv")
  des_GTEx_Colon <- des_GTEx[which(des_GTEx$group=="Colon"),]
  des_GTEx_Colon$group <- gsub("Colon","negative",des_GTEx_Colon$group)
  des_GTEx_Colon$batch <- "GTEx"
  colnames(des_GTEx_Colon) <- c("samples","group","batch")
  des_colorectal <- rbind(des_TCGA_colorectal,des_GTEx_Colon)
  unmatched <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/GTEx_unmatched_samples.txt",sep = "\t",header = TRUE)
  des_colorectal <- des_colorectal[!des_colorectal$sample %in% unmatched$unmatch_samples,]
  
  mat_raw_TCGA_GTEx_test <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/mat_TPM_TCGA_GTEx.txt",sep = "\t",header = TRUE,check.names = FALSE)
  output_res <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/colorectal_DE.txt"
  mat_raw <- mat_raw_TCGA_GTEx_test
  des <- des_colorectal
  
  result <- Wilcox_test(mat_raw,des,output_res)
  result$ensembl_gene_id <- rownames(result)
}

#Stomach cancer vs. normal tissue
{
  des_STAD <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/TCGA/STAD/TCGA_STAD_PrimaryTumor_vs_TissueNormal.csv")
  des_STAD$batch <- "TCGA"
  des_GTEx <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/GTEx/GTEx_noBlood_sample_info.csv")
  des_GTEx_Stomach <- des_GTEx[which(des_GTEx$group=="Stomach"),]
  des_GTEx_Stomach$group <- gsub("Stomach","negative",des_GTEx_Stomach$group)
  des_GTEx_Stomach$batch <- "GTEx"
  colnames(des_GTEx_Stomach) <- c("samples","group","batch")
  des_stomach <- rbind(des_STAD,des_GTEx_Stomach)
  unmatched <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/GTEx_unmatched_samples.txt",sep = "\t",header = TRUE)
  des_stomach <- des_stomach[!des_stomach$sample %in% unmatched$unmatch_samples,]
  #mat_raw_TCGA_GTEx_test <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/mat_TPM_TCGA_GTEx.txt",sep = "\t",header = TRUE,check.names = FALSE)
  output_res <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/stomach_DE.txt"
  mat_raw <- mat_raw_TCGA_GTEx_test
  des <- des_stomach
  result_stomach <- Wilcox_test(mat_raw,des,output_res)
}


#gene effect
{
  library(biomaRt)
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  gene_effect <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/CCLE/CRISPR_gene_effect.csv",check.names = FALSE,row.names = 1)
  CCLE_colon_info <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/CCLE/Colon_cell_info.csv")
  CCLE_stomach_info <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/CCLE/Stomach_cell_info.csv")
  colon_gene_effect <- as.data.frame(t(gene_effect[CCLE_colon_info[CCLE_colon_info$sample_collection_site=="Colon",]$DepMap_ID,]))
  #colon_gene_effect[is.na(colon_gene_effect)] <- 0
  stomach_gene_effect <- as.data.frame(t(gene_effect[CCLE_stomach_info[CCLE_stomach_info$sample_collection_site=="stomach",]$DepMap_ID,]))
  #stomach_gene_effect[is.na(stomach_gene_effect)] <- 0
  
  colon_gene_effect$Mean_effect <- rowMeans(colon_gene_effect,na.rm = TRUE)
  colon_gene_effect$gene_symbol <- as.character(lapply(strsplit(rownames(colon_gene_effect)," (",fixed = TRUE),function(x) x[1]))
  
  gene_names <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol","gene_biotype"),
                      filters = "hgnc_symbol",
                      values=colon_gene_effect$gene_symbol, mart= mart,useCache = FALSE)
  
  colon_gene_effect <- left_join(gene_names,colon_gene_effect,by=c("hgnc_symbol"="gene_symbol"))
  colon_FC_gene_effect <- left_join(result,colon_gene_effect,by=c("ensembl_gene_id"="ensembl_gene_id"))
  colon_FC_gene_effect[is.na(colon_FC_gene_effect)] <- 0
  
  stomach_gene_effect$Mean_effect <- rowMeans(stomach_gene_effect,na.rm = TRUE)
  stomach_gene_effect$gene_symbol <- as.character(lapply(strsplit(rownames(stomach_gene_effect)," (",fixed = TRUE),function(x) x[1]))
  
  gene_names <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol","gene_biotype"),
                      filters = "hgnc_symbol",
                      values=stomach_gene_effect$gene_symbol, mart= mart,useCache = FALSE)
  
  stomach_gene_effect <- left_join(gene_names,stomach_gene_effect,by=c("hgnc_symbol"="gene_symbol"))
  stomach_FC_gene_effect <- left_join(result,stomach_gene_effect,by=c("ensembl_gene_id"="ensembl_gene_id"))
  stomach_FC_gene_effect[is.na(stomach_FC_gene_effect)] <- 0
  #nrow(colon_FC_gene_effect[(colon_FC_gene_effect$log2FoldChange>0 & colon_FC_gene_effect$Mean_effect < 0),])
}

#plasma gene
{
  plasma_CRC <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/output/Expression_CRCvsNC_edger_exact.txt", check.names = FALSE, header = TRUE)
  plasma_STAD <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/output/Expression_STADvsNC_edger_exact.txt",check.names = FALSE, header = TRUE)
  
  plasma_CRC_upregulated <- plasma_CRC[plasma_CRC$log2FoldChange > 0,]
  plasma_CRC_upregulated$ensembl_gene_id <- as.character(lapply(strsplit(rownames(plasma_CRC_upregulated),".",fixed = TRUE),function(x) x[1]))
  plasma_STAD_upregulated <- plasma_STAD[plasma_STAD$log2FoldChange > 0,]
  plasma_STAD_upregulated$ensembl_gene_id <- as.character(lapply(strsplit(rownames(plasma_STAD_upregulated),".",fixed = TRUE),function(x) x[1]))
  
  colon_FC_gene_effect_plasmaFC <- left_join(colon_FC_gene_effect,plasma_CRC_upregulated,by=c("ensembl_gene_id"="ensembl_gene_id"))
  colon_FC_gene_effect_plasmaFC[is.na(colon_FC_gene_effect_plasmaFC)] <- 0
  stomach_FC_gene_effect_plasmaFC <- left_join(stomach_FC_gene_effect,plasma_STAD_upregulated,by=c("ensembl_gene_id"="ensembl_gene_id"))
  stomach_FC_gene_effect_plasmaFC[is.na(stomach_FC_gene_effect_plasmaFC)] <- 0
}

coupregulated_genes <- stomach_FC_gene_effect_plasmaFC[stomach_FC_gene_effect_plasmaFC$log2FoldChange.x > 0 & stomach_FC_gene_effect_plasmaFC$log2FoldChange.y > 0 & stomach_FC_gene_effect_plasmaFC$Mean_effect < 0,]$ensembl_gene_id
forenrich <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                        filters = "ensembl_gene_id",
                        values=coupregulated_genes, mart= mart,useCache = FALSE)
KEGG_stomach_res <- enrichKEGG(
  forenrich$entrezgene_id,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  #universe=as.character(background),
  minGSSize = 0,
  maxGSSize = 500,
  qvalueCutoff = 1,
  use_internal_data = FALSE)

coupregulated_genes <- colon_FC_gene_effect_plasmaFC[colon_FC_gene_effect_plasmaFC$log2FoldChange.x > 0 & colon_FC_gene_effect_plasmaFC$log2FoldChange.y > 0 & colon_FC_gene_effect_plasmaFC$Mean_effect < 0,]$ensembl_gene_id
forenrich <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id","hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values=coupregulated_genes, mart= mart,useCache = FALSE)
KEGG_colon_res <- enrichKEGG(
  forenrich$entrezgene_id,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  #universe=as.character(background),
  minGSSize = 0,
  maxGSSize = 500,
  qvalueCutoff = 1,
  use_internal_data = FALSE)

interested_genes_inpathway <- KEGG_stomach_res@result["hsa05226",]$geneID

interested_genes_inpathway <- as.character(unlist(strsplit(interested_genes_inpathway,"/",fixed = TRUE)))

getBM(attributes=c("entrezgene_id","hgnc_symbol","ensembl_gene_id"),
      filters = "entrezgene_id",
      values=interested_genes_inpathway, mart= mart,useCache = FALSE)

STAD_colon_coupregulated_genes <- intersect(stomach_FC_gene_effect_plasmaFC[stomach_FC_gene_effect_plasmaFC$log2FoldChange.x > 0 & stomach_FC_gene_effect_plasmaFC$log2FoldChange.y > 0 & stomach_FC_gene_effect_plasmaFC$Mean_effect < 0,]$ensembl_gene_id,colon_FC_gene_effect_plasmaFC[colon_FC_gene_effect_plasmaFC$log2FoldChange.x > 0 & colon_FC_gene_effect_plasmaFC$log2FoldChange.y > 0 & colon_FC_gene_effect_plasmaFC$Mean_effect < 0,]$ensembl_gene_id)
forenrich <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id","hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values=STAD_colon_coupregulated_genes, mart= mart,useCache = FALSE)
KEGG_res <- enrichKEGG(
  forenrich$entrezgene_id,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  #universe=as.character(background),
  minGSSize = 0,
  maxGSSize = 500,
  qvalueCutoff = 1,
  use_internal_data = FALSE)

genelist <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol","description","chromosome_name","start_position","end_position","strand"),
                   filters = "ensembl_gene_id",
                   values=STAD_colon_coupregulated_genes, mart= mart,useCache = FALSE)
write.csv(genelist,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/STAD_colon_coupregulated_genes.csv")
