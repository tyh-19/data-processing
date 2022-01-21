#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='Outlier detection')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
    help='input count matrix. Rows are features. Columns are samples. Tab seperated.')
parser$add_argument('-o', '--outdir', type='character', required=TRUE,
    help='output directory(full path end with /). Make sure directory exists.')

args <- parser$parse_args()

suppressPackageStartupMessages(library(progress))

#load args
  message('read count matrix: ', args$matrix)
  transfer_all <- read.csv(args$matrix,sep="\t",header = TRUE,row.names=1,stringsAsFactors = FALSE)
  gene_list <- rownames(transfer_all)
  message('feature number: ', length(gene_list))
  outdir <- args$outdir
  message('output to: ', outdir)

#load progress
pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = length(gene_list), clear = FALSE, width= 60)

#initiation
  i=1
  CRC_final <- as.data.frame(matrix(numeric(0),ncol=17))
  colnames(CRC_final) <- c("gene","FDR_up0.05","FDR_up0.01","FDR_down0.05","FDR_down0.01",
                           "CRC_up0.05","CRC_up0.05_ratio","CRC_up0.05_sample",
                           "CRC_up0.01","CRC_up0.01_ratio","CRC_up0.01_sample",
                           "CRC_down0.05","CRC_down0.05_ratio","CRC_down0.05_sample",
                           "CRC_down0.01","CRC_down0.01_ratio","CRC_down0.01_sample")
  STAD_final <- as.data.frame(matrix(numeric(0),ncol=13))
  colnames(STAD_final) <- c("gene","STAD_up0.05","STAD_up0.05_ratio","STAD_up0.05_sample",
                            "STAD_up0.01","STAD_up0.01_ratio","STAD_up0.01_sample",
                            "STAD_down0.05","STAD_down0.05_ratio","STAD_down0.05_sample",
                            "STAD_down0.01","STAD_down0.01_ratio","STAD_down0.01_sample")
  panCancer_final <- as.data.frame(matrix(numeric(0),ncol=13))
  colnames(panCancer_final) <- c("gene","panCancer_up0.05","panCancer_up0.05_ratio","panCancer_up0.05_sample",
                                 "panCancer_up0.01","panCancer_up0.01_ratio","panCancer_up0.01_sample",
                                 "panCancer_down0.05","panCancer_down0.05_ratio","panCancer_down0.05_sample",
                                 "panCancer_down0.01","panCancer_down0.01_ratio","panCancer_down0.01_sample")
  
  all_final <- as.data.frame(matrix(numeric(0),ncol=43))
  colnames(all_final) <- c("gene","FDR_up0.05","FDR_up0.01","FDR_down0.05","FDR_down0.01",
                           "CRC_up0.05","CRC_up0.05_ratio","CRC_up0.05_sample",
                           "CRC_up0.01","CRC_up0.01_ratio","CRC_up0.01_sample",
                           "CRC_down0.05","CRC_down0.05_ratio","CRC_down0.05_sample",
                           "CRC_down0.01","CRC_down0.01_ratio","CRC_down0.01_sample",
                           "gene","STAD_up0.05","STAD_up0.05_ratio","STAD_up0.05_sample",
                           "STAD_up0.01","STAD_up0.01_ratio","STAD_up0.01_sample",
                           "STAD_down0.05","STAD_down0.05_ratio","STAD_down0.05_sample",
                           "STAD_down0.01","STAD_down0.01_ratio","STAD_down0.01_sample",
                           "gene","panCancer_up0.05","panCancer_up0.05_ratio","panCancer_up0.05_sample",
                           "panCancer_up0.01","panCancer_up0.01_ratio","panCancer_up0.01_sample",
                           "panCancer_down0.05","panCancer_down0.05_ratio","panCancer_down0.05_sample",
                           "panCancer_down0.01","panCancer_down0.01_ratio","panCancer_down0.01_sample")

#FDR  
  while(i<=length(gene_list)){
  gene <- gene_list[i]
  
  transfer <- transfer_all[which(rownames(transfer_all)==gene),]

  NC <- as.data.frame(t(transfer[,grep("NC",colnames(transfer))]))
  NC$sample <- rownames(NC)
  CRC <- as.data.frame(t(transfer[,grep("CRC",colnames(transfer))]))
  CRC$sample <- rownames(CRC)
  STAD <- as.data.frame(t(transfer[,grep("STAD",colnames(transfer))]))
  STAD$sample <- rownames(STAD)
  panCancer <- as.data.frame(t(transfer[,-grep("NC",colnames(transfer))]))
  panCancer$sample <- rownames(panCancer)
  
  NC$sample <- rownames(NC)
  NC_sort <- NC[order(NC[1],decreasing = TRUE),]
  FDR_up0.05 <- NC_sort[ceiling(nrow(NC)*0.05),1]
  FDR_up0.01 <- NC_sort[ceiling(nrow(NC)*0.01),1]
  FDR_down0.05 <- NC_sort[ceiling(nrow(NC)*0.95),1]
  FDR_down0.01 <- NC_sort[ceiling(nrow(NC)*0.99),1] #ceiling is more strict for down
  
  #CRC
  {
  CRC_up0.05 <- nrow(CRC[CRC[1]>FDR_up0.05,])
  CRC_up0.05_sample <- paste(CRC[CRC[1]>FDR_up0.05,]$sample,collapse = "\n")
  CRC_up0.05_ratio <- CRC_up0.05/nrow(CRC)
  
  CRC_up0.01 <- nrow(CRC[CRC[1]>FDR_up0.01,])
  CRC_up0.01_sample <- paste(CRC[CRC[1]>FDR_up0.01,]$sample,collapse = "\n")
  CRC_up0.01_ratio <- CRC_up0.01/nrow(CRC)
  
  CRC_down0.05 <- nrow(CRC[CRC[1]<FDR_down0.05,])
  CRC_down0.05_sample <- paste(CRC[CRC[1]<FDR_down0.05,]$sample,collapse = "\n")
  CRC_down0.05_ratio <- CRC_down0.05/nrow(CRC)
  
  CRC_down0.01 <- nrow(CRC[CRC[1]<FDR_down0.01,])
  CRC_down0.01_sample <- paste(CRC[CRC[1]<FDR_down0.01,]$sample,collapse = "\n")
  CRC_down0.01_ratio <- CRC_down0.01/nrow(CRC)
  
  CRC_tmp <- data.frame(gene,FDR_up0.05,FDR_up0.01,FDR_down0.05,FDR_down0.01,
                        CRC_up0.05,CRC_up0.05_ratio,CRC_up0.05_sample,
                        CRC_up0.01,CRC_up0.01_ratio,CRC_up0.01_sample,
                        CRC_down0.05,CRC_down0.05_ratio,CRC_down0.05_sample,
                        CRC_down0.01,CRC_down0.01_ratio,CRC_down0.01_sample)
  }
  
  #STAD
  {
  STAD_up0.05 <- nrow(STAD[STAD[1]>FDR_up0.05,])
  STAD_up0.05_sample <- paste(STAD[STAD[1]>FDR_up0.05,]$sample,collapse = "\n")
  STAD_up0.05_ratio <- STAD_up0.05/nrow(STAD)
  
  STAD_up0.01 <- nrow(STAD[STAD[1]>FDR_up0.01,])
  STAD_up0.01_sample <- paste(STAD[STAD[1]>FDR_up0.01,]$sample,collapse = "\n")
  STAD_up0.01_ratio <- STAD_up0.01/nrow(STAD)
  
  STAD_down0.05 <- nrow(STAD[STAD[1]<FDR_down0.05,])
  STAD_down0.05_sample <- paste(STAD[STAD[1]<FDR_down0.05,]$sample,collapse = "\n")
  STAD_down0.05_ratio <- STAD_down0.05/nrow(STAD)
  
  STAD_down0.01 <- nrow(STAD[STAD[1]<FDR_down0.01,])
  STAD_down0.01_sample <- paste(STAD[STAD[1]<FDR_down0.01,]$sample,collapse = "\n")
  STAD_down0.01_ratio <- STAD_down0.01/nrow(STAD)
  
  STAD_tmp <- data.frame(gene,STAD_up0.05,STAD_up0.05_ratio,STAD_up0.05_sample,
                         STAD_up0.01,STAD_up0.01_ratio,STAD_up0.01_sample,
                         STAD_down0.05,STAD_down0.05_ratio,STAD_down0.05_sample,
                         STAD_down0.01,STAD_down0.01_ratio,STAD_down0.01_sample)
  }
  
  #panCancer
  {
  panCancer_up0.05 <- nrow(panCancer[panCancer[1]>FDR_up0.05,])
  panCancer_up0.05_sample <- paste(panCancer[panCancer[1]>FDR_up0.05,]$sample,collapse = "\n")
  panCancer_up0.05_ratio <- panCancer_up0.05/nrow(panCancer)
  
  panCancer_up0.01 <- nrow(panCancer[panCancer[1]>FDR_up0.01,])
  panCancer_up0.01_sample <- paste(panCancer[panCancer[1]>FDR_up0.01,]$sample,collapse = "\n")
  panCancer_up0.01_ratio <- panCancer_up0.01/nrow(panCancer)
  
  panCancer_down0.05 <- nrow(panCancer[panCancer[1]<FDR_down0.05,])
  panCancer_down0.05_sample <- paste(panCancer[panCancer[1]<FDR_down0.05,]$sample,collapse = "\n")
  panCancer_down0.05_ratio <- panCancer_down0.05/nrow(panCancer)
  
  panCancer_down0.01 <- nrow(panCancer[panCancer[1]<FDR_down0.01,])
  panCancer_down0.01_sample <- paste(panCancer[panCancer[1]<FDR_down0.01,]$sample,collapse = "\n")
  panCancer_down0.01_ratio <- panCancer_down0.01/nrow(panCancer)
  
  panCancer_tmp <- data.frame(gene,panCancer_up0.05,panCancer_up0.05_ratio,panCancer_up0.05_sample,
                              panCancer_up0.01,panCancer_up0.01_ratio,panCancer_up0.01_sample,
                              panCancer_down0.05,panCancer_down0.05_ratio,panCancer_down0.05_sample,
                              panCancer_down0.01,panCancer_down0.01_ratio,panCancer_down0.01_sample)
  }
  
  if(i%%3000==0){
    CRC_final <- rbind(CRC_final,CRC_tmp)
    STAD_final <- rbind(STAD_final,STAD_tmp)
    panCancer_final <- rbind(panCancer_final,panCancer_tmp)
    final <- cbind(CRC_final,STAD_final,panCancer_final)
    all_final <- rbind(all_final,final)
    write.csv(CRC_final,paste0(outdir,i/3000,"_CRC3000.csv"))
    write.csv(STAD_final,paste0(outdir,i/3000,"_STAD3000.csv"))
    write.csv(panCancer_final,paste0(outdir,i/3000,"_panCancer3000.csv"))
    write.csv(final,paste0(outdir,i/3000,"_3000.csv"))
    
    CRC_final <- as.data.frame(matrix(numeric(0),ncol=17))
    colnames(CRC_final) <- c("gene","FDR_up0.05","FDR_up0.01","FDR_down0.05","FDR_down0.01",
                             "CRC_up0.05","CRC_up0.05_ratio","CRC_up0.05_sample",
                             "CRC_up0.01","CRC_up0.01_ratio","CRC_up0.01_sample",
                             "CRC_down0.05","CRC_down0.05_ratio","CRC_down0.05_sample",
                             "CRC_down0.01","CRC_down0.01_ratio","CRC_down0.01_sample")
    STAD_final <- as.data.frame(matrix(numeric(0),ncol=13))
    colnames(STAD_final) <- c("gene","STAD_up0.05","STAD_up0.05_ratio","STAD_up0.05_sample",
                              "STAD_up0.01","STAD_up0.01_ratio","STAD_up0.01_sample",
                              "STAD_down0.05","STAD_down0.05_ratio","STAD_down0.05_sample",
                              "STAD_down0.01","STAD_down0.01_ratio","STAD_down0.01_sample")
    panCancer_final <- as.data.frame(matrix(numeric(0),ncol=13))
    colnames(panCancer_final) <- c("gene","panCancer_up0.05","panCancer_up0.05_ratio","panCancer_up0.05_sample",
                                   "panCancer_up0.01","panCancer_up0.01_ratio","panCancer_up0.01_sample",
                                   "panCancer_down0.05","panCancer_down0.05_ratio","panCancer_down0.05_sample",
                                   "panCancer_down0.01","panCancer_down0.01_ratio","panCancer_down0.01_sample")
    i=i+1
  } else if(i==length(gene_list)){
    CRC_final <- rbind(CRC_final,CRC_tmp)
    STAD_final <- rbind(STAD_final,STAD_tmp)
    panCancer_final <- rbind(panCancer_final,panCancer_tmp)
    final <- cbind(CRC_final,STAD_final,panCancer_final)
    all_final <- rbind(all_final,final)
    
    write.csv(CRC_final,paste0(outdir,i%/%3000+1,"_CRC",i%%3000,".csv"))
    write.csv(STAD_final,paste0(outdir,i%/%3000+1,"_STAD",i%%3000,".csv"))
    write.csv(panCancer_final,paste0(outdir,i%/%3000+1,"_panCancer",i%%3000,".csv"))
    write.csv(final,paste0(outdir,i%/%3000+1,"_",i%%3000,".csv"))
    
    i=i+1
  } else {
    CRC_final <- rbind(CRC_final,CRC_tmp)
    STAD_final <- rbind(STAD_final,STAD_tmp)
    panCancer_final <- rbind(panCancer_final,panCancer_tmp)
    i=i+1
  }
  pb$tick()
  Sys.sleep(1 / 100)
  }
  write.csv(all_final,paste0(outdir,"all.csv"))

#old version
######################################
#  pb <- progress_bar$new(
#    format = "  Processing [:bar] :percent eta: :eta",
#    total = length(gene_list), clear = FALSE, width= 60)
#  
#  i=1
#  CRC_final <- as.data.frame(matrix(numeric(0),ncol=13))
#  colnames(CRC_final) <- c("gene","FDR_up0.05","FDR_up0.01","FDR_down0.05","FDR_down0.01","CRC_up0.05","CRC_up0.05_ratio","CRC_up0.01","CRC_up0.01_ratio","CRC_down0.05","CRC_down0.05_ratio","CRC_down0.01","CRC_down0.01_ratio")
#  STAD_final <- as.data.frame(matrix(numeric(0),ncol=9))
#  colnames(STAD_final) <- c("gene","STAD_up0.05","STAD_up0.05_ratio","STAD_up0.01","STAD_up0.01_ratio","STAD_down0.05","STAD_down0.05_ratio","STAD_down0.01","STAD_down0.01_ratio")
#  panCancer_final <- as.data.frame(matrix(numeric(0),ncol=9))
#  colnames(panCancer_final) <- c("gene","panCancer_up0.05","panCancer_up0.05_ratio","panCancer_up0.01","panCancer_up0.01_ratio","panCancer_down0.05","panCancer_down0.05_ratio","panCancer_down0.01","panCancer_down0.01_ratio")
  
#  all_final <- as.data.frame(matrix(numeric(0),ncol=31))
#  colnames(all_final) <- c("gene","FDR_up0.05","FDR_up0.01","FDR_down0.05","FDR_down0.01","CRC_up0.05","CRC_up0.05_ratio","CRC_up0.01","CRC_up0.01_ratio","CRC_down0.05","CRC_down0.05_ratio","CRC_down0.01","CRC_down0.01_ratio","gene","STAD_up0.05","STAD_up0.05_ratio","STAD_up0.01","STAD_up0.01_ratio","STAD_down0.05","STAD_down0.05_ratio","STAD_down0.01","STAD_down0.01_ratio","gene","panCancer_up0.05","panCancer_up0.05_ratio","panCancer_up0.01","panCancer_up0.01_ratio","panCancer_down0.05","panCancer_down0.05_ratio","panCancer_down0.01","panCancer_down0.01_ratio")
  
#  while(i<=length(gene_list)){
#  gene <- gene_list[i]
  
#  transfer <- transfer_all[which(rownames(transfer_all)==gene),]

#  NC <- as.data.frame(t(transfer[,grep("NC",colnames(transfer))]))
#  NC$sample <- rownames(NC)
#  CRC <- as.data.frame(t(transfer[,grep("CRC",colnames(transfer))]))
#  CRC$sample <- rownames(CRC)
#  STAD <- as.data.frame(t(transfer[,grep("STAD",colnames(transfer))]))
#  STAD$sample <- rownames(STAD)
#  panCancer <- as.data.frame(t(transfer[,-grep("NC",colnames(transfer))]))
#  panCancer$sample <- rownames(panCancer)
  
#  NC$sample <- rownames(NC)
#  NC_sort <- NC[order(NC[1],decreasing = TRUE),]
#  FDR_up0.05 <- NC_sort[ceiling(nrow(NC)*0.05),1]
#  FDR_up0.01 <- NC_sort[ceiling(nrow(NC)*0.01),1]
#  FDR_down0.05 <- NC_sort[ceiling(nrow(NC)*0.95),1]
#  FDR_down0.01 <- NC_sort[ceiling(nrow(NC)*0.99),1] #ceiling is more strict for down
  
#  #CRC
#  {
#  CRC_up0.05 <- nrow(CRC[CRC[1]>FDR_up0.05,])
#  CRC_up0.05_sample <- CRC[CRC[1]>FDR_up0.05,]$sample
#  CRC_up0.05_ratio <- CRC_up0.05/nrow(CRC)
  
#  CRC_up0.01 <- nrow(CRC[CRC[1]>FDR_up0.01,])
#  CRC_up0.01_sample <- CRC[CRC[1]>FDR_up0.01,]$sample
#  CRC_up0.01_ratio <- CRC_up0.01/nrow(CRC)
  
#  CRC_down0.05 <- nrow(CRC[CRC[1]<FDR_down0.05,])
#  CRC_down0.05_sample <- CRC[CRC[1]<FDR_down0.05,]$sample
#  CRC_down0.05_ratio <- CRC_down0.05/nrow(CRC)
  
#  CRC_down0.01 <- nrow(CRC[CRC[1]<FDR_down0.01,])
#  CRC_down0.01_sample <- CRC[CRC[1]<FDR_down0.01,]$sample
#  CRC_down0.01_ratio <- CRC_down0.01/nrow(CRC)
  
#  CRC_tmp <- data.frame(gene,FDR_up0.05,FDR_up0.01,FDR_down0.05,FDR_down0.01,CRC_up0.05,CRC_up0.05_ratio,CRC_up0.01,CRC_up0.01_ratio,CRC_down0.05,CRC_down0.05_ratio,CRC_down0.01,CRC_down0.01_ratio)
#  }
  
#  #STAD
#  {
#  STAD_up0.05 <- nrow(STAD[STAD[1]>FDR_up0.05,])
#  STAD_up0.05_sample <- STAD[STAD[1]>FDR_up0.05,]$sample
#  STAD_up0.05_ratio <- STAD_up0.05/nrow(STAD)
  
#  STAD_up0.01 <- nrow(STAD[STAD[1]>FDR_up0.01,])
#  STAD_up0.01_sample <- STAD[STAD[1]>FDR_up0.01,]$sample
#  STAD_up0.01_ratio <- STAD_up0.01/nrow(STAD)
  
#  STAD_down0.05 <- nrow(STAD[STAD[1]<FDR_down0.05,])
#  STAD_down0.05_sample <- STAD[STAD[1]<FDR_down0.05,]$sample
#  STAD_down0.05_ratio <- STAD_down0.05/nrow(STAD)
  
#  STAD_down0.01 <- nrow(STAD[STAD[1]<FDR_down0.01,])
#  STAD_down0.01_sample <- STAD[STAD[1]<FDR_down0.01,]$sample
#  STAD_down0.01_ratio <- STAD_down0.01/nrow(STAD)
  
#  STAD_tmp <- data.frame(gene,STAD_up0.05,STAD_up0.05_ratio,STAD_up0.01,STAD_up0.01_ratio,STAD_down0.05,STAD_down0.05_ratio,STAD_down0.01,STAD_down0.01_ratio)
#  }
  
#  #panCancer
#  {
#  panCancer_up0.05 <- nrow(panCancer[panCancer[1]>FDR_up0.05,])
#  panCancer_up0.05_sample <- panCancer[panCancer[1]>FDR_up0.05,]$sample
#  panCancer_up0.05_ratio <- panCancer_up0.05/nrow(panCancer)
  
#  panCancer_up0.01 <- nrow(panCancer[panCancer[1]>FDR_up0.01,])
#  panCancer_up0.01_sample <- panCancer[panCancer[1]>FDR_up0.01,]$sample
#  panCancer_up0.01_ratio <- panCancer_up0.01/nrow(panCancer)
  
#  panCancer_down0.05 <- nrow(panCancer[panCancer[1]<FDR_down0.05,])
#  panCancer_down0.05_sample <- panCancer[panCancer[1]<FDR_down0.05,]$sample
#  panCancer_down0.05_ratio <- panCancer_down0.05/nrow(panCancer)
  
#  panCancer_down0.01 <- nrow(panCancer[panCancer[1]<FDR_down0.01,])
#  panCancer_down0.01_sample <- panCancer[panCancer[1]<FDR_down0.01,]$sample
#  panCancer_down0.01_ratio <- panCancer_down0.01/nrow(panCancer)
  
#  panCancer_tmp <- data.frame(gene,panCancer_up0.05,panCancer_up0.05_ratio,panCancer_up0.01,panCancer_up0.01_ratio,panCancer_down0.05,panCancer_down0.05_ratio,panCancer_down0.01,panCancer_down0.01_ratio)
#  }
  
#  if(i%%3000==0){
#    CRC_final <- rbind(CRC_final,CRC_tmp)
#    STAD_final <- rbind(STAD_final,STAD_tmp)
#    panCancer_final <- rbind(panCancer_final,panCancer_tmp)
#    final <- cbind(CRC_final,STAD_final,panCancer_final)
#    all_final <- rbind(all_final,final)
#    write.csv(CRC_final,paste0(outdir,i/3000,"_CRC3000.csv"),quote = FALSE)
#    write.csv(STAD_final,paste0(outdir,i/3000,"_STAD3000.csv"),quote = FALSE)
#    write.csv(panCancer_final,paste0(outdir,i/3000,"_panCancer3000.csv"),quote = FALSE)
#    write.csv(final,paste0(outdir,i/3000,"_3000.csv"),quote = FALSE)
    
#    CRC_final <- as.data.frame(matrix(numeric(0),ncol=13))
#    colnames(CRC_final) <- c("gene","FDR_up0.05","FDR_up0.01","FDR_down0.05","FDR_down0.01","CRC_up0.05","CRC_up0.05_ratio","CRC_up0.01","CRC_up0.01_ratio","CRC_down0.05","CRC_down0.05_ratio","CRC_down0.01","CRC_down0.01_ratio")
#    STAD_final <- as.data.frame(matrix(numeric(0),ncol=9))
#    colnames(STAD_final) <- c("gene","STAD_up0.05","STAD_up0.05_ratio","STAD_up0.01","STAD_up0.01_ratio","STAD_down0.05","STAD_down0.05_ratio","STAD_down0.01","STAD_down0.01_ratio")
#    panCancer_final <- as.data.frame(matrix(numeric(0),ncol=9))
#    colnames(panCancer_final) <- c("gene","panCancer_up0.05","panCancer_up0.05_ratio","panCancer_up0.01","panCancer_up0.01_ratio","panCancer_down0.05","panCancer_down0.05_ratio","panCancer_down0.01","panCancer_down0.01_ratio")
#    i=i+1
#  } else if(i==length(gene_list)){
#    CRC_final <- rbind(CRC_final,CRC_tmp)
#    STAD_final <- rbind(STAD_final,STAD_tmp)
#    panCancer_final <- rbind(panCancer_final,panCancer_tmp)
#    final <- cbind(CRC_final,STAD_final,panCancer_final)
#    all_final <- rbind(all_final,final)
    
#    write.csv(CRC_final,paste0(outdir,i%/%3000+1,"_CRC",i%%3000,".csv"),quote = FALSE)
#    write.csv(STAD_final,paste0(outdir,i%/%3000+1,"_STAD",i%%3000,".csv"),quote = FALSE)
#    write.csv(panCancer_final,paste0(outdir,i%/%3000+1,"_panCancer",i%%3000,".csv"),quote = FALSE)
#    write.csv(final,paste0(outdir,i%/%3000+1,"_",i%%3000,".csv"),quote = FALSE)
    
#    i=i+1
#  } else {
#    CRC_final <- rbind(CRC_final,CRC_tmp)
#    STAD_final <- rbind(STAD_final,STAD_tmp)
#    panCancer_final <- rbind(panCancer_final,panCancer_tmp)
#    i=i+1
#  }
#  pb$tick()
#  Sys.sleep(1/10)
#  }
#  write.csv(all_final,paste0(outdir,"all.csv"),quote = FALSE)

  
