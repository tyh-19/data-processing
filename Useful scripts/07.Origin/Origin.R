##deconvolution by cibersortx (previous)
{
{
  #install.packages("remotes")
  #remotes::install_github("icbi-lab/immunedeconv")
  library(ggplot2)
  library(immunedeconv)
  library(tidyverse)
}
##GTEx
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/GTEx/TSS/")
  
  signature <- read.csv("signature_matrix_TSS0.9_TPM_top20.txt",sep="\t",header = TRUE,row.names = 1)
  col_annotation <- read.csv("sample_ids_annotation.txt",sep = "\t", header = FALSE,row.names = 1)
  row_annotation <- read.csv("signature_genes_forplot.txt", sep = "\t", header = FALSE, row.names = 1)
  colnames(row_annotation) <- c("Tissue_specific_gene")
  colnames(col_annotation) <- c("Tissue_type","Tissue_type_specific")
  rownames(col_annotation)<- gsub("-",".",rownames(col_annotation))
  
  #signature matrix plot
  bk = unique(c(seq(-0.5,0.5, length=100)))
  pheatmap(
    signature,
    breaks = bk,
    annotation_col = col_annotation,
    annotation_row = row_annotation,
    scale = "row",
    cluster_cols = FALSE,cluster_rows = FALSE,
    show_colnames=FALSE, 
    show_rownames=FALSE,
    #cluster_cols = FALSE,
    colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 5)
  
  signature_matrix <- read.csv("signature_matrix_TSS0.9_TPM_top20_median.txt",sep = "\t", header = T,row.names = 1)
  pheatmap(
    signature_matrix,
    breaks = bk,
    #annotation_col = col_annotation,
    annotation_row = row_annotation,
    scale = "row",
    cluster_cols = FALSE,cluster_rows = FALSE,
    show_colnames=TRUE, 
    show_rownames=FALSE,
    #cluster_cols = FALSE,
    colorRampPalette(c("blue","white","red"))(100),
    fontsize_col = 10)
  
  
  #preprocess for cibersortx
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/GTEx/forcibersortx_TSS/")
  signature <- read.csv("signature_genes_forplot.txt",sep="\t",header = FALSE,row.names = 1)
  mixture <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/pico_tissue_TPM.txt",sep="\t",header = T, row.names = 1)
  
  #rownames(mixture) <- as.character(lapply(strsplit(rownames(mixture),".",fixed = TRUE),function(x) x[1]))
  rownames(signature) <- as.character(lapply(strsplit(rownames(signature),".",fixed = TRUE),function(x) x[1]))
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  annotations <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), 
                       filters = "ensembl_gene_id", 
                       values=rownames(signature), mart= mart,useCache = FALSE)
  
  
  
  #get RNAs
  j=1
  pathway_gene_count={}
  while(j<=nrow(annotations)){
    target <- annotations[j,2]
    gene_symbol <- annotations[j,1]
    if(length(grep(target,rownames(mixture)))==0) {
      #print(paste0("No ",target," in this dataset."))
      temp <- as.data.frame(array(,dim=c(1,ncol(mixture))))
      temp[1,] <- 0
      rownames(temp) <- target
      colnames(temp) <- colnames(mixture)
    } else {
      #temp <- mixture[which(rownames(mixture)==target),]
      #temp <- mixture[grep(target,rownames(mixture),fixed=TRUE),]  #for ensg
      temp <- mixture[grep(target,rownames(mixture),fixed=TRUE),] #for gene symbol
      rownames(temp) <- target
    }
    pathway_gene_count <- rbind(pathway_gene_count,temp)
    j=j+1
  }
  
  #pathway_gene_count$symbol <- annotations$hgnc_symbol
  write.csv(pathway_gene_count,"pico_tissue_forGTEx_signature_TSS.csv")
  #write.csv(signature,"signature_matrix_DE_top50_TPM_rmdup_ENSG.csv",quote = F)
}

#Cibersort_locally
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/GTEx/forcibersortx/")
  source('/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/GTEx/forcibersortx/Cibersort.R')  
  res_cibersort <- CIBERSORT("signature_matrix_DE_top50_TPM_rmdup_ENSG_median.txt","Plasma_forGTEx_signature.txt",perm = 1000, QN = F)
  write.table(res_cibersort, file="res_cibersort.txt", sep="\t", col.names=T, row.names=F, quote=F)
  
  library(ggplot2)
  library(tidyverse)
  colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  
  my_theme <- function(){
    theme(panel.grid = element_blank(),       # 网格线
          panel.border = element_blank(),     # 面板边框
          legend.position="right",            # legend位置
          legend.text = element_text(size=8), # legend内容大小
          legend.title = element_text(size=8),# legend标题大小
          axis.line = element_line(size=1),   # 坐标轴线
          text = element_text(family="Times"),# 文本字体
          axis.text.y = element_text(size = 8,face='bold',color='black'),# y轴标签样式
          axis.text.x = element_text(size = 8,face='bold',color='black',angle=90,hjust=1),        # x轴标签样式，angle=45 倾斜 45 度
          axis.title = element_text(size=10,face="bold"),  # 轴标题
          plot.title = element_text(hjust=0.5,size=10))    # 距，设置绘图区域距离边的据类，上、右、下、左
  }  
  
  p1 <- res_cibersort[,1:22] %>% reshape2::melt() %>%
    ggplot(aes(x=Var1,y=value,fill=Var2)) +
    geom_bar(stat='identity') +
    coord_flip()  +
    scale_fill_manual(values =colour ) +
    theme_bw()+ theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+my_theme()
  
  pdf("cibersort.pdf")
  p1
}

#stackplot
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/LM22/LM22_deconvolution/PBMC/")
  composition <- read.csv("tissue_origin_PBMC.csv",header = TRUE)
  
  stackplot <- melt(composition)
  
  #boxplot
  {
    stackplot$Type <- as.character(lapply(strsplit(as.character(stackplot$variable),".",fixed = TRUE),function(x) x[3])) 
    stackplot$Type <- gsub("T","Tumor Tissue", stackplot$Type)
    stackplot$Type <- gsub("N","Normal Tissue", stackplot$Type)
    
    my_comparisons <- list(c("Normal Tissue","Tumor Tissue"))
    p <- ggplot(stackplot[which(stackplot$Average=="Blood-Vessel"),],aes(x=Type,y=value,fill=Type))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
      scale_fill_brewer(palette="Blues") +
      theme_bw()+
      theme(#legend.position="right",
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))
    
    m=-1
    if(m>0){
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                method.args = list(alternative = "greater"),
                                label = "p.signif"
      )+labs(x="",y="Fraction (%)",title="wilcox.test.greater", face="bold",fill="Type")
      p
    } else if(m==0){
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                method.args = list(alternative = "two.sided"),
                                label = "p.signif"
      )+labs(x="",y="Fraction (%)",title="wilcox.test.twoside", face="bold",fill="Type")
      p
    } else {
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                method.args = list(alternative = "less"),
                                label = "p.signif"
      )+labs(x="",y="Fraction (%)",title="Blood-Vessel\nwilcox.test.less", face="bold",fill="Type")
      p
    }
  }
  
  #plot plasma
  stackplot$Average <- factor(stackplot$Average,levels = c("Blood-Whole-Blood","Blood-Vessel","Colon","Esophagus-Muscularis","Liver","Lung","Stomach"))
  stackplot$variable <- factor(stackplot$variable,levels = c("Colorectum","Stomach","Lung","Liver","Esophagus","Healthy"))
  
  #plot tissue
  comp_tmp <- stackplot[which(stackplot$Average=="T cells CD8"),]
  
  stackplot$Average <- factor(stackplot$Average,levels = unique(sort(stackplot$Average)))
  stackplot$variable <- factor(stackplot$variable,levels = comp_tmp[order(comp_tmp$value,decreasing = TRUE),2])
  
  ggplot(stackplot,aes(x=stackplot$variable,y=stackplot$value*100,fill = stackplot$Average)) + geom_bar(stat = "identity", width=0.5, col='black') +
    theme_bw()+
    theme(#legend.position="bottom",
      legend.position="right",
      panel.grid=element_blank(),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20,angle = 90,hjust = 1,vjust =0.5),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    ylab("Fraction(%)")+
    xlab("")+
    labs(fill="Tissue")+
    #geom_vline(aes(xintercept=6.5))+
    scale_y_continuous(breaks = c(0,25,50,75,100),labels = c("0","25","50","75","100"),expand = c(0,0),limits = c(0,103))
  scale_fill_nejm(alpha = 1)
}

#make median matrix and tissue specific score
{
  sample_info <- as.data.frame(list(c("a","b","c","d","e","f","g","h","i"),c("A","A","A","B","B","B","C","C","C")))
  colnames(sample_info) <- c("sample","group")
  mat <- as.data.frame(list(c(100,10,1),c(110,8,2),c(120,6,3),c(10,110,1),c(8,120,2),c(6,100,3),c(1,10,110),c(2,8,120),c(3,6,100)))
  colnames(mat) <- c("a","b","c","d","e","f","g","h","i")
  rownames(mat) <- c("gene1","gene2","gene3")
  
  group <- unique(sort(sample_info$group))
  i=1
  while(i<=length(group)){
    sample <- sample_info[which(sample_info$group==group[i]),]$sample
    median <- as.data.frame(rowMedians(as.matrix(mat[,sample])))
    colnames(median) <- group[i]
    rownames(median) <- rownames(mat)
    if(i==1){
      median_matrix <- median
    } else {
      median_matrix <- cbind(median_matrix,median)
    }
    i=i+1
  }
  
  j=1
  sum <- as.data.frame(rowSums(median_matrix))
  while(j<=ncol(median_matrix)){
    TSS <- median_matrix[,j]/sum
    colnames(TSS) <- colnames(median_matrix)[j]
    if(j==1){
      TSS_matrix <- TSS
    } else {
      TSS_matrix <- cbind(TSS_matrix,TSS)
    }
    j=j+1
  }
  
}

#LM22
#input mixture pre-process
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/")
  mixture <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/pico_tissue_TPM.txt",sep="\t",header = T, row.names = 1)
  signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/LM22.txt",sep="\t",header = T, row.names = 1)
  
  #rownames(mixture) <- as.character(lapply(strsplit(rownames(mixture),".",fixed = TRUE),function(x) x[1])) 
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  annotations <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), 
                       filters = "hgnc_symbol", 
                       values=rownames(signature), mart= mart,useCache = FALSE)
  
  #get RNAs
  j=1
  pathway_gene_count={}
  while(j<=nrow(annotations)){
    target <- annotations[j,2]
    gene_symbol <- annotations[j,1]
    if(length(grep(target,rownames(mixture)))==0) {
      #print(paste0("No ",target," in this dataset."))
      temp <- as.data.frame(array(,dim=c(1,ncol(mixture))))
      temp[1,] <- 0
      rownames(temp) <-target
      colnames(temp) <- colnames(mixture)
    } else {
      #temp <- mixture[which(rownames(mixture)==target),]
      #temp <- mixture[grep(target,rownames(mixture),fixed=TRUE),]  #for ensg
      temp <- mixture[grep(target,rownames(mixture),fixed=TRUE),] #for gene symbol
      rownames(temp) <- target
    }
    pathway_gene_count <- rbind(pathway_gene_count,temp)
    j=j+1
  }
  
  pathway_gene_count$symbol <- annotations$hgnc_symbol
  write.csv(pathway_gene_count,"Tissue_forLM22.csv")
}
}

#deconvolution
{
library(matrixStats)
library(Biobase)
library(pheatmap)
library(dplyr)
#function
mean_TSS <- function(mat,sample_info,TPM_cutoff,TSS_cutoff){
  sample_info_ex <- as.data.frame(list(c("a","b","c","d","e","f","g","h","i"),c("A","A","A","B","B","B","C","C","C")))
  colnames(sample_info_ex) <- c("sample","group")
  mat_ex <- as.data.frame(list(c(100,10,1),c(110,8,2),c(120,6,3),c(10,110,1),c(8,120,2),c(6,100,3),c(1,10,110),c(2,8,120),c(3,6,100)))
  colnames(mat_ex) <- c("a","b","c","d","e","f","g","h","i")
  rownames(mat_ex) <- c("gene1","gene2","gene3")
  
  message("Example mat:")
  print(mat_ex)
  message("Example sample_info:")
  print(sample_info_ex)
  message("Example TPM_cutoff: 20")
  message("Example TSS_cutoff: 0.95")
  
  group <- unique(sort(sample_info$group))
  i=1
  while(i<=length(group)){
    sample <- sample_info[which(sample_info$group==group[i]),]$sample
    median <- as.data.frame(rowMeans(as.matrix(mat[,sample])))
    colnames(median) <- group[i]
    rownames(median) <- rownames(mat)
    if(i==1){
      median_matrix <- median
    } else {
      median_matrix <- cbind(median_matrix,median)
    }
    i=i+1
  }

  median_matrix <- median_matrix[which(rowMins(as.matrix(median_matrix)) > TPM_cutoff),]
  
  j=1
  sum <- as.data.frame(rowSums(median_matrix))
  while(j<=ncol(median_matrix)){
    TSS <- median_matrix[,j]/sum
    colnames(TSS) <- colnames(median_matrix)[j]
    if(j==1){
      TSS_matrix <- TSS
    } else {
      TSS_matrix <- cbind(TSS_matrix,TSS)
    }
    j=j+1
  }
  TSG_plot <- TSS_matrix[which(rowMaxs(as.matrix(TSS_matrix)) > TSS_cutoff),]
  
  k=1
  TSG_plot$group <- NA 
  while(k<=length(group)){
    if(length(grep("TRUE",TSG_plot[,group[k]] > TSS_cutoff))>0){
      TSG_plot[(TSG_plot[,group[k]] > TSS_cutoff),]$group <- group[k]
      k=k+1
    } else {
      k=k+1
    }
  }
  TSG_plot <- TSG_plot[order(TSG_plot$group),]
  signature_matrix <- median_matrix[rownames(TSG_plot),]
  signature_matrix <- signature_matrix[order(TSG_plot$group),]
  rownames(signature_matrix) <- as.character(lapply(strsplit(rownames(signature_matrix),".",fixed = TRUE),function(x) x[1]))
  write.table(signature_matrix,"signature_genes.txt",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
  return(TSG_plot)
}

#blueprint: blood
{
  #summarize rawcount matrix
  {
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin")
  records <- readxl::read_xlsx("Blueprint_RNAseq.xlsx",sheet = "Sheet1")
  unique(records$`Cell type`)

  i=1
  files <- records$`count file name`
  count <- read.table(paste0("Blueprint/",files[1]),header = TRUE)
  all <- data.frame(count$gene_id)
  while(i<=length(files)){
    message(paste0(i,"/",length(files)))
    count <- read.table(paste0("Blueprint/",files[i]),header = TRUE)
    tmp <- count[,c("gene_id","TPM")]
    colnames(tmp) <- c("gene_id",files[i])
    all <- dplyr::left_join(all,tmp,by=c("count.gene_id"="gene_id"))
    i=i+1
    }      
  count_matrix <- all
  #write.table(all,"Blueprint.txt",sep = "\t",quote = FALSE,row.names = FALSE)
  
  colnames(count_matrix) <- gsub(".1..","(1).",colnames(count_matrix),fixed=TRUE)
  colnames(count_matrix) <- gsub(".2..","(2).",colnames(count_matrix),fixed=TRUE)
  colnames(count_matrix) <- gsub(".3..","(3).",colnames(count_matrix),fixed=TRUE)
  colnames(count_matrix) <- gsub(".4..","(4).",colnames(count_matrix),fixed=TRUE)
  colnames(count_matrix) <- gsub(".5..","(5).",colnames(count_matrix),fixed=TRUE)

  write.table(count_matrix[,-grep(").",colnames(count_matrix))],"Blueprint_TPM.sorted.txt",sep = "\t",quote = FALSE, row.names = TRUE)
  }
  
  #build signature matrix
  { 
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin")
    count_matrix <- read.csv("Blueprint_TPM.sorted.txt",sep = "\t",header = TRUE,row.names = 1)
    sample_information <- readxl::read_xlsx("Blueprint_RNAseq.sorted.xlsx")
    sample_info <- as.data.frame(sample_information[-grep("x",sample_information$`Category[customized]`),c("count file name","Category[customized]")])
    #sample_info <- as.data.frame(sample_information[-which(sample_information$`Cell type...12`=="-"),c("count file name","Cell type...12")])
    colnames(sample_info) <- c("sample","group")
    sample_info$sample <- gsub("-",".",sample_info$sample,fixed = TRUE)
    mat <- count_matrix[,sample_info$sample]
    TPM_cutoff <- 0
    TSS_cutoff <- 0.70
    
    TSG_plot <- mean_TSS(mat,sample_info,TPM_cutoff,TSS_cutoff)
    
    col_annotation <- sample_info
    rownames(col_annotation) <- col_annotation$sample
    col_annotation$sample <- NULL
    #col_annotation <- col_annotation[,c("group")]
    
    bk = unique(c(seq(0,1, length=100)))
    pheatmap(
      mat[rownames(TSG_plot),order(sample_info$group)],
      breaks = bk,
      annotation_col = col_annotation,
      #annotation_row = row_annotation,
      scale = "row",
      cluster_cols = FALSE,cluster_rows = FALSE,
      show_colnames=FALSE, 
      show_rownames=FALSE,
      colorRampPalette(c("blue","white","red"))(100),
      fontsize_row = 5,
      clustering_distance_cols = "euclidean")
    
    #plot signature gene matrix
    library(biomaRt)
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    
    id <- as.character(lapply(strsplit(rownames(TSG_plot),"\\."),function(x) x[1]))
    
    transversion <- data.frame("raw"=rownames(TSG_plot),"id"=id)
    
    gene_names <- getBM(attributes=c("ensembl_gene_id", "external_gene_name","gene_biotype"),
                        filters = "ensembl_gene_id",
                        values=transversion$id, mart= mart,useCache = FALSE)
    
    annotated <- left_join(transversion,gene_names, by= c("id"="ensembl_gene_id"))
    
    rownames(TSG_plot) <- paste(annotated$gene_biotype,annotated$raw,annotated$external_gene_name,sep = "|")
    
    bk = unique(c(seq(0,1, length=100)))
    pheatmap(
      TSG_plot[order(TSG_plot$group),-which(colnames(TSG_plot)=="group")],
      breaks = bk,
      #annotation_col = col_annotation,
      #annotation_row = row_annotation,
      scale = "row",
      cluster_cols = FALSE,cluster_rows = FALSE,
      show_colnames=TRUE, 
      show_rownames=FALSE,
      colorRampPalette(c("blue","white","red"))(100),
      fontsize_row = 5,
      clustering_distance_cols = "euclidean")
    }
  
  #local cibersort
  source('Cibersort.R')
  gene_expression_matrix <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/ipico_gencode_TPM.txt",sep = "\t")
  rownames(gene_expression_matrix) <- as.character(lapply(strsplit(rownames(gene_expression_matrix),"|",fixed = TRUE),function(x) x[1]))
  write.table(gene_expression_matrix,"TPM/ipico_gencode_TPM.txt",sep = "\t", quote = FALSE, row.names = TRUE,col.names = TRUE)
  
  res_cibersort <- CIBERSORT("signature/blood_signature_genes.txt","TPM/ipico_gencode_TPM.txt",perm = 1000, QN = F)
  write.table(res_cibersort, file="output/res_cibersort_ipico.txt", sep="\t", col.names=T, row.names=T, quote=F)
  
  
  library(ggplot2)
  library(tidyverse)
  #full plot
  {
  res_cibersort<- read.csv("Summarized_results.csv",header = TRUE)
  colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  
  my_theme <- function(){
    theme(panel.grid = element_blank(),       # 网格线
          panel.border = element_blank(),     # 面板边框
          legend.position="right",            # legend位置
          legend.text = element_text(size=8,face='bold',color='black'), # legend内容大小
          legend.title = element_blank(),# legend标题大小
          axis.line.x = element_line(size=0.5),   # 坐标轴线
          axis.line.y = element_blank(),
          text = element_text(family="Arial"),# 文本字体
          axis.text.y = element_text(size = 8,face='bold',color='black'),# y轴标签样式
          axis.text.x = element_text(size = 8,face='bold',color='black',angle=90,hjust=1),        # x轴标签样式，angle=45 倾斜 45 度
          axis.title = element_blank(),  # 轴标题
          plot.title = element_text(hjust=0.5,size=10))    # 距，设置绘图区域距离边的据类，上、右、下、左
  }  
  
  p1 <- res_cibersort[order(res_cibersort$sample.type),1:9] %>% reshape2::melt()
  p1$X <- factor(p1$X,level =  rev(c("GSE142987","GSE174302","multiomics","ipico_Plasma",
                                 "ipico_EV","exoRBase","GSE133684",
                                 "GSE68086","GSE89843","GSE156902",
                                 "in-housed PBMC","GSE182522",
                                 "GSE181465","in-housed tissue")))
  p1$variable <- factor(p1$variable,level =  c("B.cell","T.cell","NK.cell",
                                                   "Monocytes","Granulocytes",
                                                   "Megakaryocytes","Erythrocytes",
                                                   "Non.hematopoietic"))
  plot <- ggplot(p1,aes(x=X,y=value,fill=variable)) +
    geom_bar(stat='identity') +
    coord_flip()  +
    scale_fill_manual(values =colour ) +
    #scale_fill_jco()+
    scale_y_continuous(limits = c(0,1.01),expand = c(0,0)) + 
    theme_bw()+ theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+my_theme()
  
  ggsave(plot = plot,filename = "Composition.pdf",device = "pdf",height = 3.09,width = 4.76)
  
  }
  
  #simplified plot
  {
    res_cibersort<- read.csv("Summarized_results.Simplified.csv",header = TRUE)
    colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
    
    my_theme <- function(){
      theme(panel.grid = element_blank(),       # 网格线
            panel.border = element_blank(),     # 面板边框
            legend.position="right",            # legend位置
            legend.text = element_text(size=8,face='bold',color='black'), # legend内容大小
            legend.title = element_blank(),# legend标题大小
            axis.line.x = element_line(size=0.5),   # 坐标轴线
            axis.line.y = element_blank(),
            text = element_text(family="Arial"),# 文本字体
            axis.text.y = element_text(size = 8,face='bold',color='black'),# y轴标签样式
            axis.text.x = element_text(size = 8,face='bold',color='black',angle=90,hjust=1),        # x轴标签样式，angle=45 倾斜 45 度
            axis.title = element_blank(),  # 轴标题
            plot.title = element_text(hjust=0.5,size=10))    # 距，设置绘图区域距离边的据类，上、右、下、左
    }  
    
    p1 <- res_cibersort[order(res_cibersort$sample.type),1:6] %>% reshape2::melt()
    p1$X <- factor(p1$X,level =  rev(c("GSE142987","GSE174302","multiomics","ipico_Plasma",
                                       "ipico_EV","exoRBase","GSE133684",
                                       "GSE68086","GSE89843","GSE156902",
                                       "in-housed PBMC","GSE182522",
                                       "GSE181465","in-housed tissue")))
    p1$variable <- factor(p1$variable,level =  c("PBMC","Granulocytes",
                                                 "Megakaryocytes","Erythrocytes",
                                                 "Non.hematopoietic"))
    plot <- ggplot(p1,aes(x=X,y=value,fill=variable)) +
      geom_bar(stat='identity') +
      coord_flip()  +
      scale_fill_manual(values =colour ) +
      #scale_fill_jco()+
      scale_y_continuous(limits = c(0,1.01),expand = c(0,0)) + 
      theme_bw()+ theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+my_theme()
    ggsave(plot = plot,filename = "Composition.simplified.pdf",device = "pdf",height = 3.09,width = 4.76)
  }
}

#tissue: HPA+GETx consensus 55 tissue
{
  tissue <- as.data.frame(read_tsv("rna_tissue_consensus.tsv"))
  
  tissue$GENE <- paste(tissue$Gene,tissue$`Gene name`,sep = "|")
  
  tissue_matrix <- as.data.frame(array(dim=c(length(unique(tissue$GENE)),length(unique(tissue$Tissue)))))
  colnames(tissue_matrix) <- unique(tissue$Tissue)
  rownames(tissue_matrix) <- unique(tissue$GENE)
  
  i=1
  while(i<=nrow(tissue)){
    tissue_matrix[tissue[i,"GENE"],tissue[i,"Tissue"]] <- tissue[i,"nTPM"]
    i=i+1
  }
  
  tissue_TSS <- tissue_matrix
  tissue_TSS[1:nrow(tissue_TSS),1:ncol(tissue_TSS)] <- NA
  tissue_sum <- rowSums(tissue_matrix)
  i=1
  while(i<=nrow(tissue_matrix)){
    tissue_TSS[i,] <- tissue_matrix[i,]/tissue_sum[i]
    i=i+1
  }
  
  View(tissue_matrix[which(rowMaxs(as.matrix(tissue_TSS)) > 0.99),])
  TSG_plot <- tissue_matrix[which(rowMaxs(as.matrix(tissue_TSS)) > 0.99),]
 
  
  bk = unique(c(seq(-0.5,0.5, length=100)))
  pheatmap(
    TSG_plot,
    breaks = bk,
    #annotation_col = col_annotation,
    #annotation_row = row_annotation,
    scale = "row",
    cluster_cols = TRUE,cluster_rows = TRUE,
    show_colnames=TRUE, 
    show_rownames=FALSE,
    colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 5,
    clustering_distance_cols = "euclidean")
}

#tissue: single cell type tissue 
{
  tissue <- as.data.frame(read_tsv("rna_single_cell_type_tissue.tsv"))
  
  tissue$GENE <- paste(tissue$Gene,tissue$`Gene name`,sep = "|")
  
  tissue_matrix <- as.data.frame(array(dim=c(length(unique(tissue$GENE)),length(unique(tissue$Tissue)))))
  colnames(tissue_matrix) <- unique(tissue$Tissue)
  rownames(tissue_matrix) <- unique(tissue$GENE)
  
  i=1
  while(i<=nrow(tissue)){
    tissue_matrix[tissue[i,"GENE"],tissue[i,"Tissue"]] <- tissue[i,"nTPM"]
    i=i+1
  }
  
  
}

#Blueprint+GTEx+TCGA (need larger than 6G memory)
{
  library(matrixStats)
  library(Biobase)
  library(pheatmap)
  library(dplyr)
  #function
  mean_TSS <- function(mat,sample_info,TPM_cutoff,TSS_cutoff){
    sample_info_ex <- as.data.frame(list(c("a","b","c","d","e","f","g","h","i"),c("A","A","A","B","B","B","C","C","C")))
    colnames(sample_info_ex) <- c("sample","group")
    mat_ex <- as.data.frame(list(c(100,10,1),c(110,8,2),c(120,6,3),c(10,110,1),c(8,120,2),c(6,100,3),c(1,10,110),c(2,8,120),c(3,6,100)))
    colnames(mat_ex) <- c("a","b","c","d","e","f","g","h","i")
    rownames(mat_ex) <- c("gene1","gene2","gene3")
    
    message("Example mat:")
    print(mat_ex)
    message("Example sample_info:")
    print(sample_info_ex)
    message("Example TPM_cutoff: 20")
    message("Example TSS_cutoff: 0.95")
    
    group <- unique(sort(sample_info$group))
    i=1
    while(i<=length(group)){
      sample <- sample_info[which(sample_info$group==group[i]),]$sample
      median <- as.data.frame(rowMeans(as.matrix(mat[,sample])))
      colnames(median) <- group[i]
      rownames(median) <- rownames(mat)
      if(i==1){
        median_matrix <- median
      } else {
        median_matrix <- cbind(median_matrix,median)
      }
      i=i+1
    }
    
    median_matrix <- median_matrix[which(rowMins(as.matrix(median_matrix)) > TPM_cutoff),]
    
    j=1
    sum <- as.data.frame(rowSums(median_matrix))
    while(j<=ncol(median_matrix)){
      TSS <- median_matrix[,j]/sum
      colnames(TSS) <- colnames(median_matrix)[j]
      if(j==1){
        TSS_matrix <- TSS
      } else {
        TSS_matrix <- cbind(TSS_matrix,TSS)
      }
      j=j+1
    }
    TSG_plot <- TSS_matrix[which(rowMaxs(as.matrix(TSS_matrix)) > TSS_cutoff),]
    
    k=1
    TSG_plot$group <- NA 
    while(k<=length(group)){
      if(length(grep("TRUE",TSG_plot[,group[k]] > TSS_cutoff))>0){
        TSG_plot[(TSG_plot[,group[k]] > TSS_cutoff),]$group <- group[k]
        k=k+1
      } else {
        k=k+1
      }
    }
    TSG_plot <- TSG_plot[order(TSG_plot$group),]
    signature_matrix <- median_matrix[rownames(TSG_plot),]
    signature_matrix <- signature_matrix[order(TSG_plot$group),]
    rownames(signature_matrix) <- as.character(lapply(strsplit(rownames(signature_matrix),".",fixed = TRUE),function(x) x[1]))
    write.table(signature_matrix,"signature_genes.txt",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
    return(TSG_plot)
  }
  
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/Blood_GTEx_TCGA")
  #read in TPM
  {
  Blood <- read.csv("Blueprint_TPM.sorted.txt",sep = "\t",header = TRUE,row.names = 1)
  rownames(Blood) <- as.character(lapply(strsplit(rownames(Blood),".",fixed = TRUE),function(x) x[1]))
  
  Tissue <- read.csv("GTEx_ensembl_gene_tpm.txt",sep = "\t",header = TRUE, row.names = 1)
  rownames(Tissue) <- as.character(lapply(strsplit(rownames(Tissue),".",fixed = TRUE),function(x) x[1]))
  
  Tumor_COAD <- read.csv("TCGA_COAD_count_matrix_TPM.txt",sep = "\t", header = TRUE, row.names = 1,check.names = FALSE)
  Tumor_READ <- read.csv("TCGA_READ_count_matrix_TPM.txt",sep = "\t", header = TRUE, row.names = 1,check.names = FALSE)
  Tumor_STAD <- read.csv("TCGA_STAD_count_matrix_TPM.txt",sep = "\t", header = TRUE, row.names = 1,check.names = FALSE)
  
  rownames(Tumor_COAD) <- as.character(lapply(strsplit(rownames(Tumor_COAD),"|",fixed = TRUE),function(x) x[1]))
  rownames(Tumor_READ) <- as.character(lapply(strsplit(rownames(Tumor_READ),"|",fixed = TRUE),function(x) x[1]))
  rownames(Tumor_STAD) <- as.character(lapply(strsplit(rownames(Tumor_STAD),"|",fixed = TRUE),function(x) x[1]))
  }
  
  #combine Blood GTEx and TCGA
  {
    Blood$Blood <- rownames(Blood)
    Tissue$Tissue <- rownames(Tissue)
    Tumor_COAD$Tumor_COAD <- rownames(Tumor_COAD)
    Tumor_READ$Tumor_READ <- rownames(Tumor_READ)
    Tumor_STAD$Tumor_STAD <- rownames(Tumor_STAD)
    
    Tumor_CRC <- left_join(Tumor_COAD,Tumor_READ,by=c("Tumor_COAD"="Tumor_READ"))
    Tumor_GI <- left_join(Tumor_CRC, Tumor_STAD, by=c("Tumor_COAD"="Tumor_STAD"))
    
    Blood_tumor <- left_join(Blood,Tumor_GI,by = c("Blood"="Tumor_COAD"))
    
    Blood_tumor_tissue <- left_join(Blood_tumor,Tissue, by = c("Blood"="Tissue"))
    
    rownames(Blood_tumor_tissue) <- Blood_tumor_tissue$Blood
    Blood_tumor_tissue <- Tumor_GI[,-which(colnames(Blood_tumor_tissue)=="Blood")]
  }
  
  #get signature genes and average expression
  {
    count_matrix <- Blood_tumor_tissue
    
    sample_info <- read.csv("Blood_GTEx_Tumor_sample_info.csv")
    
    mat <- count_matrix[,sample_info$sample]
    TPM_cutoff <- 0
    TSS_cutoff <- 0.70
    
    TSG_plot <- mean_TSS(mat,sample_info,TPM_cutoff,TSS_cutoff)
  }
  }
}

#result presentation
#devtools::install_github("ricardo-bion/ggradar", 
#                         dependencies = TRUE)
library(ggplot2)
library(ggradar)
library(scales)

composition <-read.csv("output/res_cibersort_multiomics_Blood_Tissue.txt",sep="\t",header=TRUE,row.names=1)
composition<-composition[grep("pico",rownames(composition)),]
composition<-composition[-grep("mix..pico",rownames(composition)),]
composition <- composition[-which(rownames(composition)=="CRC.PKU.29.pico"),]

composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),".",fixed = TRUE),function(x) x[1]))
composition$group<-gsub("NC","HD",composition$group)

composition$group <- factor(composition$group,levels = c("HD","CRC","STAD"))
my_comparisons <- list(c("HD","CRC"),c("HD","STAD"),c("CRC","STAD"))

ggplot(composition,aes(x=composition$group,y=Stomach,fill = composition$group))+
  #ggplot(forplot,aes(x=group,y=value,fill = group))+
  geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
  geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
  #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
  #scale_fill_brewer(palette="Blues") +
  #ylim(0,25)+
  theme_bw()+
  xlab("")+
  #ylab("")+
  #ylab(colnames(plot)[1])+
  theme(#legend.position="right",
    legend.position="none",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(size=1, colour = "black"),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
    axis.text.y = element_text(face="bold",  color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "two.sided",paired = TRUE),
                     label = "p.signif",
                     size = 10,
                     vjust = 0.5)

mtcars %>%
  add_rownames( var = "group" ) %>%
  mutate_each(funs(rescale), -group) %>%
  tail(4) %>% dplyr::select(1:10) -> mtcars_radar

ggradar(mtcars_radar) 

radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=mean)
radar <- as.tibble(radar)
radar <- radar[,c(1,rev(order(radar[1,-1]))+1)]
radar <- radar[,-which(colnames(radar)=="RMSE")]
radar <- radar[,-which(colnames(radar)=="Correlation")]
radar <- radar[,-which(colnames(radar)=="P.value")]
#radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
ggradar(radar[,c("Group","COAD","READ","Colon","STAD","Stomach")],grid.min = 0,grid.mid = 0.00125, grid.max = 0.0025,
        values.radar = c("0%", "0.125%", "0.25%"),
        font.radar = "Arial",
        group.point.size = 1,
        group.line.width = 0.5)

ggradar(radar[,c("Group","B.cell","T.cell","NK.cell")],grid.min = 0,grid.mid = 0.0003125, grid.max = 0.000625,
        values.radar = c("0%", "0.03125%", "0.0625%"),
        font.radar = "Arial",
        group.point.size = 1,
        group.line.width = 0.5)

ggradar(radar[,1:6],grid.min = 0,grid.mid = 0.5, grid.max = 1,
        plot.extent.x.sf = 1.5, plot.extent.y.sf = 2,
        values.radar = c("0%", "50%", "100%"),
        font.radar = "Arial",
        group.point.size = 1,
        group.line.width = 0.5,
        group.colours = c(CRC="#FCB514",STAD="#CD9B1D",HD="#87CEEB",GIC="#FCB514"))

#LM22
{
composition <-read.csv("output/res_cibersort_multiomics_LM22.txt",sep="\t",header=TRUE,row.names=1)
composition<-composition[grep("pico",rownames(composition)),]
composition<-composition[-grep("mix..pico",rownames(composition)),]
composition<-composition[-grep("CRC.PKU.29.pico",rownames(composition)),]

composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),".",fixed = TRUE),function(x) x[1]))
composition$group<-gsub("NC","HD",composition$group)
composition$group<-gsub("CRC","GIC",composition$group)
composition$group<-gsub("STAD","GIC",composition$group)

radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=mean)
radar <- as.tibble(radar)
radar <- radar[,c(1,rev(order(radar[2,-1]))+1)]
radar <- radar[,-which(colnames(radar)=="RMSE")]
radar <- radar[,-which(colnames(radar)=="Correlation")]
radar <- radar[,-which(colnames(radar)=="P.value")]
colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
#radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
ggradar(radar,grid.min = 0,grid.mid = 0.1, grid.max = 0.2,
        values.radar = c("", "10%", "20%"),
        plot.extent.x.sf = 1.5, plot.extent.y.sf = 2,
        font.radar = "Arial",
        group.point.size = 2,
        group.line.width = 1,
        group.colours = c(CRC="#FCB514",STAD="#CD9B1D",HD="#87CEEB",GIC="#FCB514"))
}

#noCancer
{
  composition <-read.csv("output/res_cibersort_multiomics_noCancer.simplified.txt",sep="\t",header=TRUE,row.names=1)
  composition<-composition[grep("pico",rownames(composition)),]
  composition<-composition[-grep("mix..pico",rownames(composition)),]
  composition<-composition[-grep("CRC.PKU.29.pico",rownames(composition)),]
  
  composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),".",fixed = TRUE),function(x) x[1]))
  composition$group<-gsub("NC","HD",composition$group)
  composition$group<-gsub("CRC","GIC",composition$group)
  composition$group<-gsub("STAD","GIC",composition$group)
  
  radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=mean)
  radar <- as.tibble(radar)
  radar <- radar[,c(1,rev(order(radar[2,-1]))+1)]
  radar <- radar[,-which(colnames(radar)=="RMSE")]
  radar <- radar[,-which(colnames(radar)=="Correlation")]
  radar <- radar[,-which(colnames(radar)=="P.value")]
  colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
  #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
  ggradar(radar,grid.min = 0,grid.mid = 0.5, grid.max = 1,
          plot.extent.x.sf = 1.5, plot.extent.y.sf = 2,
          values.radar = c("", "50%", "100%"),
          font.radar = "Arial",
          group.point.size = 2,
          group.line.width = 1,
          group.colours = c(CRC="#FCB514",STAD="#CD9B1D",HD="#87CEEB",GIC="#FCB514"))
}
