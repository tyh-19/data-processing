#! /usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(biomaRt))
parser <- ArgumentParser(description='Signature Matrix v2')
parser$add_argument('-m', '--matrix', type='character', required=TRUE,
                    help='Input count matrix. Rows are genes(TPM). Columns are samples.')
parser$add_argument('-s', '--sample_info', type='character', required=TRUE,
                    help='Sample information for each sample. CSV file with 2 columns. Column 1: sample, Column 2: group')
parser$add_argument('-e', '--TPM_cutoff', type='double', required=FALSE,default=1,
                    help='TPM cutoff, maximum expressionn larger than this cutoff. Default: 1.')
parser$add_argument('-ts','--TSS_cutoff', type='double', required=FALSE,default=0.5,
                    help='TSS cutoff, minimum tissue specific score larger than this cutoff. Default: 0.5.')
parser$add_argument('-n','--Number', type='integer', required=FALSE,default=20,
                    help='Signature gene number. Default: 20.')
parser$add_argument('-o', '--outdir', type='character', required=TRUE,
                    help='output directory')
args <- parser$parse_args()

#cl <- makeCluster(1) #not to overload your computer
#registerDoParallel(cl)

message('read count matrix: ', args$matrix)
message('output: ', args$outdir)
message('Start read in files at: ',Sys.time())
#function
mean_TSS <- function(mat,sample_info,TPM_cutoff=1,TSS_cutoff=0.5,Num=20,outdir){
  sample_info_ex <- as.data.frame(list(c("a","b","c","d","e","f","g","h","i"),c("A","A","A","B","B","B","C","C","C")))
  colnames(sample_info_ex) <- c("sample","group")
  mat_ex <- as.data.frame(list(c(100,10,1),c(110,8,2),c(120,6,3),c(10,110,1),c(8,120,2),c(6,100,3),c(1,10,110),c(2,8,120),c(3,6,100)))
  colnames(mat_ex) <- c("a","b","c","d","e","f","g","h","i")
  rownames(mat_ex) <- c("gene1","gene2","gene3")
  
  message("Example mat:")
  print(mat_ex)
  message("Example sample_info:")
  print(sample_info_ex)
  message("Example TPM_cutoff: 1")
  message("Example TSS_cutoff: 0.5")
  message("Example Number: 20")
  
  group <- unique(sort(sample_info$group))
  i=1
  while(i<=length(group)){
    sample <- sample_info[which(sample_info$group==group[i]),]$sample
    median <- as.data.frame(rowMeans(as.matrix(mat[,as.character(sample)])))
    colnames(median) <- group[i]
    rownames(median) <- rownames(mat)
    if(i==1){
      median_matrix <- median
    } else {
      median_matrix <- cbind(median_matrix,median)
    }
    i=i+1
  }
  
  median_matrix <- median_matrix[which(rowMaxs(as.matrix(median_matrix)) > TPM_cutoff),]
  
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

  n=1
  TSG_plot={}
  while(n<=ncol(TSS_matrix)){
    test_sorted <- TSS_matrix[order(TSS_matrix[,n],decreasing = TRUE),]
    test_cutoff <- test_sorted[test_sorted[,n]>0.5,]
    TSG_plot <- rbind(TSG_plot,test_cutoff[1:Num,])
    n=n+1
  }  

  TSG_plot <- TSG_plot[which(rowMaxs(as.matrix(TSG_plot)) > TSS_cutoff),]
  
  k=1
  TSG_plot$group <- NA 
  while(k<=length(group)){
    if(length(grep("TRUE",TSG_plot[,as.character(group[k])] > TSS_cutoff))>0){
      TSG_plot[(TSG_plot[,as.character(group[k])] > TSS_cutoff),]$group <- as.character(group[k])
      k=k+1
    } else {
      k=k+1
    }
  }
  
  #TSG_plot <- TSG_plot[order(TSG_plot$group),]
  signature_matrix <- median_matrix[rownames(TSG_plot),]
  signature_matrix <- signature_matrix[order(TSG_plot$group),]
  rownames(signature_matrix) <- as.character(lapply(strsplit(rownames(signature_matrix),".",fixed = TRUE),function(x) x[1]))
  write.table(signature_matrix,paste0(outdir,"/","signature_genes_",TPM_cutoff,"_",TSS_cutoff,"_",Num,".txt"),quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
  row_annotation <- data.frame(rownames(TSG_plot),TSG_plot[,which(colnames(TSG_plot)=="group")])
  colnames(row_annotation) <- c("gene","group")
  rownames(row_annotation) <- row_annotation$gene
  row_annotation$gene <- NULL
  write.table(row_annotation,paste0(outdir,"/","gene_list_",TPM_cutoff,"_",TSS_cutoff,"_",Num,".txt"),quote=FALSE,sep="\t",row.names = TRUE, col.names = TRUE)
  return(TSG_plot)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

outdir <- args$outdir
count_matrix <- read.table(args$matrix,sep = "\t",header = TRUE, row.names = 1, check.names = FALSE)
sample_info <- read.csv(args$sample_info,header = TRUE)

all_sample <- as.data.frame(colnames(count_matrix))
colnames(all_sample) <- "sample"
all_sample$exist <- "Y"
sample_info <- left_join(sample_info,all_sample, by=c("sample"="sample"))
sample_info <- sample_info[which(sample_info$exist=="Y"),]
sample_info$exist <- NULL

mat <- count_matrix[,as.character(sample_info$sample)]
TPM_cutoff <- args$TPM_cutoff
TSS_cutoff <- args$TSS_cutoff
Num <- args$Number

message('TPM cutoff: ', TPM_cutoff)
message('TSS cutoff: ', TSS_cutoff)
message('Signature gene number: ', Num)
message('Start calculating tissue specific genes at: ',Sys.time())

TSG_plot <- mean_TSS(mat,sample_info,TPM_cutoff,TSS_cutoff,Num,outdir)

message('End calculating tissue specific genes at: ',Sys.time())
message('Start ploting signature genes at: ',Sys.time())

col_annotation <- sample_info
rownames(col_annotation) <- col_annotation$sample
col_annotation <- col_annotation[order(sample_info$group),]
col_annotation$sample <- NULL
#col_annotation <- col_annotation[,c("group")]

mat_plot <- mat[rownames(TSG_plot),order(sample_info$group)]
bk = unique(c(seq(0,1, length=100)))
plot_sample <- pheatmap(
  mat_plot[order(TSG_plot$group),],
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
save_pheatmap_pdf(plot_sample, paste0(outdir,"/","signature_genes_samples_",TPM_cutoff,"_",TSS_cutoff,"_",Num,".pdf"))
dev.off()

#plot signature gene matrix
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#id <- as.character(lapply(strsplit(rownames(TSG_plot),"\\."),function(x) x[1]))
#transversion <- data.frame("raw"=rownames(TSG_plot),"id"=id)
#gene_names <- getBM(attributes=c("ensembl_gene_id", "external_gene_name","gene_biotype"),
#                    filters = "ensembl_gene_id",
#                    values=transversion$id, mart= mart,useCache = FALSE)
#annotated <- left_join(transversion,gene_names, by= c("id"="ensembl_gene_id"))
#rownames(TSG_plot) <- paste(annotated$gene_biotype,annotated$raw,annotated$external_gene_name,sep = "|")

row_annotation <- data.frame(rownames(TSG_plot),TSG_plot[,which(colnames(TSG_plot)=="group")])
colnames(row_annotation) <- c("gene","group")
rownames(row_annotation) <- row_annotation$gene
row_annotation$gene <- NULL
#row_annotation <- row_annotation[order(TSG_plot$group),]
#head(row_annotation[order(TSG_plot$group),])

bk = unique(c(seq(0,1, length=100)))
plot_tissue <- pheatmap(
  TSG_plot[order(TSG_plot$group),-which(colnames(TSG_plot)=="group")],
  breaks = bk,
  #annotation_col = col_annotation,
  annotation_row = row_annotation,
  scale = "row",
  cluster_cols = FALSE,cluster_rows = FALSE,
  show_colnames=TRUE, 
  show_rownames=FALSE,
  colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 5,
  clustering_distance_cols = "euclidean")

save_pheatmap_pdf(plot_tissue, paste0(outdir,"/","signature_genes_tissues_",TPM_cutoff,"_",TSS_cutoff,"_",Num,".pdf"))
dev.off()
