#function
{
  #function for enrichment analysis for RNA expression: clusterprofiler
  KEGG_GO <- function(genes,output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC){
    ensembl_id <- as.character(lapply(strsplit(genes,"\\|"),function(x) x[1]))
    forenrich <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                            filters = "ensembl_gene_id",
                            values=ensembl_id, mart= mart,useCache = FALSE)
    if(length(which(forenrich$entrezgene_id==""))==0){
      forenrich <- forenrich
    } else {
      forenrich <- forenrich[-which(forenrich$entrezgene_id==""),]
    }
    forenrich <- forenrich$entrezgene_id
    forenrich <- as.character(na.omit(forenrich))
    
    #KEGG
    {
      KEGG_res <- enrichKEGG(
        forenrich,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1,
        use_internal_data = FALSE)
      KEGG_output <- KEGG_res@result
      write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
    }
    #GO_BP
    {
      library(clusterProfiler)
      GO_res <- enrichGO(
        forenrich,
        'org.Hs.eg.db',
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output <- GO_res@result
      write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
    }
    #GO_CC
    {
      library(clusterProfiler)
      GO_res <- enrichGO(
        forenrich,
        'org.Hs.eg.db',
        ont = "CC",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output <- GO_res@result
      write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
    }
    #GO_MF
    {
      library(clusterProfiler)
      GO_res <- enrichGO(
        forenrich,
        'org.Hs.eg.db',
        ont = "MF",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output <- GO_res@result
      write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
    }
  }
}


##PCAWG, expression/splicing
{
  #expression
  {
    outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/FDR/" #expression
    expression_all <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/plasma/TPM_intron_spanning_noMTRNA_PassedQC.txt",sep = "\t",header = TRUE, row.names = 1)
    
    transfer_all <- expression_all[1:10,]
    gene_list <- rownames(transfer_all)
  }
  
  #splicing
  {
    outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/FDR/" #splicing
    Splicing_A3SS <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/A3SS_JC_inc_level.txt",sep="\t",header = TRUE,row.names=1)
    Splicing_A5SS <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/A5SS_JC_inc_level.txt",sep="\t",header = TRUE,row.names=1)
    Splicing_MXE <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/MXE_JC_inc_level.txt",sep="\t",header = TRUE,row.names=1)
    Splicing_SE <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/SE_JC_inc_level.txt",sep="\t",header = TRUE,row.names=1)
    Splicing_RI <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/RI_JC_inc_level.txt",sep="\t",header = TRUE,row.names=1)
    
    splicing <- rbind(Splicing_A3SS,Splicing_A5SS,Splicing_MXE,Splicing_SE,Splicing_RI)
    
    transfer_all <- splicing
    
    #get all passed FDR<0.05 splicing events
    Splicing_A3SS_stats <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/A3SS_JC_stats.txt",sep="\t",header = TRUE,row.names=1)
    Splicing_A5SS_stats <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/A5SS_JC_stats.txt",sep="\t",header = TRUE,row.names=1)
    Splicing_MXE_stats <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/MXE_JC_stats.txt",sep="\t",header = TRUE,row.names=1)
    Splicing_SE_stats <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/SE_JC_stats.txt",sep="\t",header = TRUE,row.names=1)
    Splicing_RI_stats <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/RI_JC_stats.txt",sep="\t",header = TRUE,row.names=1)
    
    splicing_stats <- rbind(Splicing_A3SS_stats,Splicing_A5SS_stats,Splicing_MXE_stats,Splicing_SE_stats,Splicing_RI_stats)
    
    gene_list <- rownames(splicing_stats[splicing_stats$FDR<0.05,])
    
    Splicing_FDR0.05 <- {}
    i=1
    while (i<=length(gene_list)){
      print(i)
      tmp <- splicing[which(rownames(splicing)==gene_list[i]),]
      Splicing_FDR0.05 <- rbind(Splicing_FDR0.05,tmp)
      i=i+1
    }
    colnames(Splicing_FDR0.05) <- colnames(splicing)
    Splicing_FDR0.05[is.na(Splicing_FDR0.05)] <- 0
    colnames(Splicing_FDR0.05) <- gsub("X.data.taoyuhuan.projects.exOmics_RNA.multiomics_paired.output.Intron.spanning.","",fixed = TRUE,colnames(Splicing_FDR0.05))
    colnames(Splicing_FDR0.05) <- gsub(".intron.spanning.bam","",fixed = TRUE,colnames(Splicing_FDR0.05))
    colnames(Splicing_FDR0.05) <- gsub(".","-",fixed = TRUE,colnames(Splicing_FDR0.05))
    write.table(Splicing_FDR0.05,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/splicing/Splicing_FDR.txt",quote = FALSE,sep="\t")
  }
  
  #APA
  {
    outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/APA/FDR/" #APA
    APA <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/APA/PDUI.txt",sep="\t",header = TRUE,row.names=1)
    transfer_all <- APA
    gene_list <- rownames(APA)
  }
  
  #TE
  {
    outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/TE/FDR/" #TE
    TE <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/TE/multiomics_paired_TE_TPM_PassedQC.txt",sep="\t",header = TRUE,row.names=1,stringsAsFactors = FALSE)
    TE <- TE[-1,]
    transfer_all <- as.data.frame(lapply(TE,as.numeric))
    rownames(transfer_all) <- rownames(TE)
    gene_list <- rownames(transfer_all)
  }
  
  #chimeric
  {
    outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Chimeric RNA/FDR/" #chimeric
    chimeric <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Chimeric RNA/multiomics_paired_chimeric_PassedQC.txt",sep="\t",header = TRUE,row.names=1,stringsAsFactors = FALSE)
    transfer_all <- chimeric
    gene_list <- rownames(transfer_all)
  }
  
  #Editing
  {
    outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/REDI/FDR/" #Editing
    editing <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/REDI/editing_all_PassedQC.txt",sep="\t",header = TRUE,row.names=1,stringsAsFactors = FALSE)
    editing[is.na(editing)] <- 0
    transfer_all <- editing
    gene_list <- rownames(editing)
  }
  
  #ASE #too large, must run on cluster
  {
    outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Paired/ASE/FDR" #ASE
    ASE_CRC_NC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Paired/ASE/CRC_NC.txt",sep="\t",header = TRUE,row.names=1,stringsAsFactors = FALSE)
    ASE_STAD_NC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Paired/ASE/STAD_NC.txt",sep="\t",header = TRUE,row.names=1,stringsAsFactors = FALSE)
    
    ASE_STAD <- ASE_STAD_NC[,grep("STAD",colnames(ASE_STAD_NC))]
    
    ASE_CRC_NC$index <- rownames(ASE_CRC_NC)
    ASE_STAD$index <- rownames(ASE_STAD)
    ASE_CRC_STAD_NC <- full_join(ASE_CRC_NC,ASE_STAD,by=c("index"="index"))
    
    grep("index",colnames(ASE_CRC_STAD_NC),value=TRUE)
    rownames(ASE_CRC_STAD_NC) <- ASE_CRC_STAD_NC$index
    ASE_CRC_STAD_NC <- ASE_CRC_STAD_NC[,-grep("index",colnames(ASE_CRC_STAD_NC))]
    ASE_CRC_STAD_NC[is.na(ASE_CRC_STAD_NC)] <- 0
    colnames(ASE_CRC_STAD_NC) <- gsub(".","-",fixed = TRUE,colnames(ASE_CRC_STAD_NC))
    write.table(ASE_CRC_STAD_NC,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Paired/ASE/CRC_STAD_NC.txt",quote = FALSE,sep="\t")
    
    transfer_all <- ASE_CRC_STAD_NC
    gene_list <- rownames(ASE_CRC_STAD_NC)
  }
  
  #SNP #too large, must run on cluster
  {
    outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Paired/SNP/FDR" #SNP
    SNP_CRC_NC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Paired/SNP/CRC_NC.txt",sep="\t",header = TRUE,row.names=1,stringsAsFactors = FALSE)
    SNP_STAD_NC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Paired/SNP/STAD_NC.txt",sep="\t",header = TRUE,row.names=1,stringsAsFactors = FALSE)
    
    SNP_STAD <- SNP_STAD_NC[,grep("STAD",colnames(SNP_STAD_NC))]
    SNP_CRC_NC$index <- rownames(SNP_CRC_NC)
    SNP_STAD$index <- rownames(SNP_STAD)
    SNP_CRC_STAD_NC <- full_join(SNP_CRC_NC,SNP_STAD,by=c("index"="index"))
    grep("index",colnames(SNP_CRC_STAD_NC),value=TRUE)
    rownames(SNP_CRC_STAD_NC) <- SNP_CRC_STAD_NC$index
    SNP_CRC_STAD_NC <- SNP_CRC_STAD_NC[,-grep("index",colnames(SNP_CRC_STAD_NC))]
    SNP_CRC_STAD_NC[is.na(SNP_CRC_STAD_NC)] <- 0
    colnames(SNP_CRC_STAD_NC) <- gsub(".","-",fixed = TRUE,colnames(SNP_CRC_STAD_NC))
    write.table(SNP_CRC_STAD_NC,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Paired/SNP/CRC_STAD_NC.txt",quote = FALSE,sep="\t")
    
    transfer_all <- chimeric
    gene_list <- rownames(transfer_all)
  }
  
  #Alt.promoter #too large, need to run on cluster
  {
    outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Alt.promoter/FDR" #Alt.promoter
    Alt.promoter <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Alt.promoter/TPM-by-tx_multiomics_paired.txt",sep="\t",header = TRUE,row.names=1,stringsAsFactors = FALSE)
    transfer_all <- Alt.promoter
    gene_list <- rownames(transfer_all)
  }
  
  #microbe
  {
    outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/microbe/FDR/" #microbe
    microbe_G <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/microbe/microbe_G_level.txt",sep="\t",header = TRUE,row.names=1)
    microbe_G[is.na(microbe_G)] <- 0
    transfer_all <- microbe_G
    gene_list <- rownames(microbe_G)
    
  }
  
  #Outlier Finder
  {
    pb <- progress_bar$new(
      format = "  Processing [:bar] :percent eta: :eta",
      total = length(gene_list), clear = FALSE, width= 60)
    
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
  }
  
  #PCAWG plot
  {
    #load functions
    {
      #select best performance
      select_bestperformance <- function(data,level,group){
        print("data: rows are chr_position, or chr|position, chromsome and position are seperated by '_' or '|'. Or gene ensembl ids.")
        print("level: confident level 0.05 or 0.01.")
        print("group: CRC or STAD or panCancer")
        
        data$BP_level <- NA
        data$BP_level_ratio <- NA
        data$BP_level_sample <- NA
        
        i=1
        while(i<=nrow(data)){
          up_ratio <- data[i,which(colnames(data)==paste0(group,"_up",level,"_ratio"))]
          down_ratio <- data[i,which(colnames(data)==paste0(group,"_down",level,"_ratio"))]
          up_n <- which(colnames(data)==paste0(group,"_up",level,"_sample"))
          down_n <- which(colnames(data)==paste0(group,"_down",level,"_sample"))
          up_sample <- as.character(data[i,as.numeric(up_n)])
          down_sample  <- as.character(data[i,as.numeric(down_n)])
          up <- data[i,which(colnames(data)==paste0(group,"_up",level))]
          down <- data[i,which(colnames(data)==paste0(group,"_down",level))]
          
          if(up_ratio>=down_ratio){
            data[i,which(colnames(data)=="BP_level")] <- up
            data[i,which(colnames(data)=="BP_level_ratio")] <- up_ratio
            data[i,which(colnames(data)=="BP_level_sample")] <- up_sample
          } else {
            data[i,which(colnames(data)=="BP_level")] <- down
            data[i,which(colnames(data)=="BP_level_ratio")] <- down_ratio
            data[i,which(colnames(data)=="BP_level_sample")] <- down_sample
          }
          
          #show remain jobs
          if(i%%10000==0){
            print(paste0("Finished: ",round(i/nrow(data)*100,2),"%"))
          } else if(i==nrow(data)) {
            print("Finished!")
          } 
          i=i+1
          
        }
        
        cutoff_tmp <- paste0("BP_",group,"_",level)
        colnames(data) <- c(colnames(data[,1:(ncol(data)-3)]),cutoff_tmp,paste0(cutoff_tmp,"_ratio"),paste0(cutoff_tmp,"_sample"))
        return(data)
      }
      
      #chromosome location2ensembl
      location2ensembl <- function(data,candidate_list){
        print("data: rows are chr_position, or chr|position, chromsome and position are seperated by '_' or '|'.")
        print("list must contain 5 columns: Gene, ensembl, chromosome, start, end.")
        
        #prepare input matrix rownames
        rownames(data) <- gsub("|","_",fixed=TRUE,rownames(data))
        data$gene <- as.character(data$gene)
        data$chromosome <- as.character(lapply(strsplit(as.character(rownames(data)),"_",fixed = TRUE),function(x) x[1]))
        data$position <- as.integer(lapply(strsplit(as.character(rownames(data)),"_",fixed = TRUE),function(x) x[2]))
        
        i=1
        while(i<=nrow(candidate_list)){
          candidate <- candidate_list$ensembl[i]
          name <- candidate_list$Gene[i]
          chromosome <- as.character(candidate_list$chromosome[i])
          start <- candidate_list$start[i]
          end <- candidate_list$end[i]
          
          if(nrow(data[data$chromosome %in% chromosome & data$position>=start & data$position<=end & !is.na(data$position),])==0) {
            print(paste0(i,": No ",candidate,":",name," in this dataset."))
            #data[data$chromosome %in% chromosome & data$position>=start & data$position<=end & !is.na(data$position), which(colnames(data)=="gene")] <- paste(name,"|",candidate)
            i=i+1
          } else {
            data[data$chromosome %in% chromosome & data$position>=start & data$position<=end & !is.na(data$position), which(colnames(data)=="gene")] <- paste(name,"|",candidate)
            i=i+1
          } 
        }
        data <- data[grep("ENSG",data$gene),]
        return(data)
      }
      
      #subset by ensembl id
      subset_row <- function(data,candidate_list,cutoff){
        print("data: rows are gene ensembl ids.")
        print("list must contain 2 columns: Gene, ensembl.")
        i=1
        output_matrix <- {}
        while(i<=nrow(candidate_list)){
          candidate <- candidate_list$ensembl[i]
          name <- candidate_list$Gene[i]
          if(length(grep(candidate,data$gene))==0) {
            #print(paste0("No ",candidate," in this dataset."))
            temp <- as.data.frame(array(,dim=c(1,ncol(data))))
            temp[1,] <- NA
            colnames(temp) <- colnames(data)
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          } else if(length(grep(candidate,data$gene))==1) {
            temp <- data[grep(candidate,data$gene,fixed=TRUE),]
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          } else {
            ##max
            minidata <- data[grep(candidate,data$gene),]
            minidata <- minidata[order(minidata[,cutoff],decreasing = TRUE),]
            temp <- minidata[1,]
            #average
            #
            #union
            #
            #take 1st as represent
            #temp <- data[grep(candidate,data$gene,fixed=TRUE)[1],]
            
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          }
        }
        return(output_matrix)
      }
      
      #chimeric subset by gene_symbol
      chimeric_subset_row <- function(data,candidate_list,cutoff){
        print("data: rows are gene enembl ids.")
        print("list must contain 2 columns: Gene, ensembl.")
        i=1
        output_matrix <- {}
        while(i<=nrow(candidate_list)){
          candidate <- candidate_list$ensembl[i]
          name <- candidate_list$Gene[i]
          
          toMatch <- c(paste0("^",name,"\\|"),paste0("\\|",name,"\\|"))
          pattern <- paste(toMatch,collapse="|")
          
          if(length(grep(pattern,data$gene))==0) {
            #print(paste0("No ",candidate," in this dataset."))
            temp <- as.data.frame(array(,dim=c(1,ncol(data))))
            temp[1,] <- NA
            colnames(temp) <- colnames(data)
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          } else if(length(grep(pattern,data$gene))==1) {
            temp <- data[grep(pattern,data$gene),]
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          } else {
            ##max
            minidata <- data[grep(pattern,data$gene),]
            minidata <- minidata[order(minidata[,cutoff],decreasing = TRUE),]
            temp <- minidata[1,]
            ##average
            #
            ##union
            #
            ##take 1st as represent
            #temp <- data[grep(pattern,data$gene)[1],]
            
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          }
        }
        return(output_matrix)
      }
    }
    
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/PCAWG_20211028/")
    candidate_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Colorectal Cancer Genes.csv",header = TRUE)
    #choose a confidence level
    level <- 0.01 #c(0.05,0.01)
    #choose an interested cancer type
    group <- "panCancer" #c("CRC","STAD","panCancer")
    #choose a cutoff, then subset
    cutoff <- "BP_panCancer_0.01" #c("BP_CRC_0.05","BP_CRC_0.01","BP_STAD_0.05","BP_STAD_0.01","BP_panCancer_0.05","BP_panCancer_0.01","CRC_up0.05","CRC_up0.01","CRC_down0.05","CRC_down0.01","STAD_up0.05","STAD_up0.01","STAD_down0.05","STAD_down0.01","panCancer_up0.05","panCancer_up0.01","panCancer_down0.05","panCancer_down0.01")
    
    #give the total sample number
    
    #CRC sample number
    #{
    #sampleN_RNA <- 23
    #sampleN_DNA <- 23
    #}
    #STAD sample number
    #{
    #  sampleN_RNA <- 30
    #  sampleN_DNA <- 30
    #}
    
    #panCancer sample number
    {
      sampleN_RNA <- 53
      sampleN_DNA <- 53
    }
    
    #read in feature outliers
    {
      #Alt.promoter
      Alt <- read.csv("Altpromoter/all.csv",header = TRUE)
      Alt <- Alt[,-1]
      rownames(Alt) <- Alt$gene
      Alt <- select_bestperformance(Alt,level,group)
      Alt_targetgenes <- subset_row(Alt,candidate_list,cutoff)
      
      #Expression
      Expression <- read.csv("Expression/all.csv",header = TRUE)
      Expression <- Expression[,-1]
      rownames(Expression) <- Expression$gene
      Expression <- select_bestperformance(Expression,level,group)
      Expression_targetgenes <- subset_row(Expression,candidate_list,cutoff)
      
      #Splicing
      Splicing <- read.csv("Splicing/all.csv",header = TRUE)
      Splicing <- Splicing[,-1]
      rownames(Splicing) <- Splicing$gene
      Splicing <- select_bestperformance(Splicing,level,group)
      Splicing_targetgenes <- subset_row(Splicing,candidate_list,cutoff)
      
      #APA
      #APA <- read.csv("APA/all.csv",header = TRUE)
      #APA <- APA[,-1]
      #rownames(APA) <- APA$gene
      #APA <- select_bestperformance(APA,level,group)
      #APA_targetgenes <- subset_row(APA,candidate_list,cutoff)
      
      #chimeric
      chimeric <- read.csv("chimeric/all.csv",header = TRUE)
      chimeric <- chimeric[,-1]
      rownames(chimeric) <- chimeric$gene
      chimeric <- select_bestperformance(chimeric,level,group)
      chimeric_targetgenes <- chimeric_subset_row(chimeric,candidate_list,cutoff)
      
      #Editing
      Editing <- read.csv("Editing/all.csv",header = TRUE)
      Editing <- Editing[,-1]
      rownames(Editing) <- Editing$gene
      Editing <- location2ensembl(Editing,candidate_list)
      Editing <- select_bestperformance(Editing,level,group)
      Editing_targetgenes <- subset_row(Editing,candidate_list,cutoff)
      
      #ASE
      ASE <- read.csv("ASE/all.csv",header = TRUE)
      ASE <- ASE[,-1]
      rownames(ASE) <- ASE$gene
      ASE <- location2ensembl(ASE,candidate_list)
      ASE <- select_bestperformance(ASE,level,group)
      ASE_targetgenes <- subset_row(ASE,candidate_list,cutoff)
      
      #SNP
      SNP <- read.csv("SNP/all.csv",header = TRUE)
      SNP <- SNP[,-1]
      rownames(SNP) <- SNP$gene
      SNP <- location2ensembl(SNP,candidate_list)
      SNP <- select_bestperformance(SNP,level,group)
      SNP_targetgenes <- subset_row(SNP,candidate_list,cutoff)
      
      #DNA-CNV
      DNA_CNV <- read.csv("DNA-CNV/all.csv",header = TRUE)
      DNA_CNV <- DNA_CNV[,-1]
      rownames(DNA_CNV) <- DNA_CNV$gene
      DNA_CNV <- select_bestperformance(DNA_CNV,level,group)
      DNA_CNV_targetgenes <- subset_row(DNA_CNV,candidate_list,cutoff)
      
      #DNA-nucleosome
      DNA_nucleosome <- read.csv("DNA-nucleosome/all.csv",header = TRUE)
      DNA_nucleosome <- DNA_nucleosome[,-1]
      rownames(DNA_nucleosome) <- DNA_nucleosome$gene
      DNA_nucleosome <- select_bestperformance(DNA_nucleosome,level,group)
      DNA_nucleosome_targetgenes <- subset_row(DNA_nucleosome,candidate_list,cutoff)
      
      #MeDIP
      DNA_MeDIP <- read.csv("MeDIP/all.csv",header = TRUE)
      DNA_MeDIP <- DNA_MeDIP[,-1]
      rownames(DNA_MeDIP) <- DNA_MeDIP$gene
      DNA_MeDIP <- select_bestperformance(DNA_MeDIP,level,group)
      DNA_MeDIP_targetgenes <- subset_row(DNA_MeDIP,candidate_list,cutoff)
      
    }
    
    #plot
    {
      #summarize outliers (RNA+DNA merged barplot for complementary)
      {
        #input features
        Alt_cutoff <- Alt_targetgenes[grep(cutoff,colnames(Alt_targetgenes))]
        Alt_cutoff[,3] <- gsub(".pico", "", Alt_cutoff[,3])
        
        Expression_cutoff <- Expression_targetgenes[grep(cutoff,colnames(Expression_targetgenes))]
        Expression_cutoff[,3] <- gsub(".pico", "", Expression_cutoff[,3])
        
        Splicing_cutoff <- Splicing_targetgenes[grep(cutoff,colnames(Splicing_targetgenes))]
        Splicing_cutoff[,3] <- gsub(".pico", "", Splicing_cutoff[,3])
        
        #APA_cutoff <- APA_targetgenes[grep(cutoff,colnames(APA_targetgenes))]
        #APA_cutoff[,3] <- gsub(".pico", "", APA_cutoff[,3])
        
        chimeric_cutoff <- chimeric_targetgenes[grep(cutoff,colnames(chimeric_targetgenes))]
        chimeric_cutoff[,3] <- gsub(".pico", "", chimeric_cutoff[,3])
        
        Editing_cutoff <- Editing_targetgenes[grep(cutoff,colnames(Editing_targetgenes))]
        Editing_cutoff[,3] <- gsub(".pico", "", Editing_cutoff[,3])
        
        ASE_cutoff <- ASE_targetgenes[grep(cutoff,colnames(ASE_targetgenes))]
        ASE_cutoff[,3] <- gsub(".pico", "", ASE_cutoff[,3])
        
        SNP_cutoff <- SNP_targetgenes[grep(cutoff,colnames(SNP_targetgenes))]
        SNP_cutoff[,3] <- gsub(".pico", "", SNP_cutoff[,3])
        
        DNA_CNV_cutoff <- DNA_CNV_targetgenes[grep(cutoff,colnames(DNA_CNV_targetgenes))]
        DNA_CNV_cutoff[,3] <- gsub(".wgs", "", DNA_CNV_cutoff[,3])
        
        DNA_nucleosome_cutoff <- DNA_nucleosome_targetgenes[grep(cutoff,colnames(DNA_nucleosome_targetgenes))]
        DNA_nucleosome_cutoff[,3] <- gsub(".wgs", "", DNA_nucleosome_cutoff[,3])
        
        DNA_MeDIP_cutoff <- DNA_MeDIP_targetgenes[grep(cutoff,colnames(DNA_MeDIP_targetgenes))]
        DNA_MeDIP_cutoff[,3] <- gsub(".me", "", DNA_MeDIP_cutoff[,3])
        
        features <- cbind(Alt_cutoff,Expression_cutoff,Splicing_cutoff,
                          chimeric_cutoff,Editing_cutoff,ASE_cutoff,SNP_cutoff,
                          DNA_CNV_cutoff,DNA_nucleosome_cutoff,DNA_MeDIP_cutoff) ## this sequential is important! It determines the label seqeuntial mannually set at melt step!
        
        features$UnionN <- NA 
        
        #Union number of input features
        i=1
        while(i<=nrow(features)){
          x <- features[i,grep(paste0(cutoff,"_sample"),colnames(features))]
          print(x)
          y <- lapply(x,as.character)
          z <- paste(y,sep = "\n",collapse = "\n")
          u <- unique(unlist(strsplit(z,"\n")))
          u <- u[u!=""]
          u <- u[u!="NA"]
          features[i,which(colnames(features)=="UnionN")] <- length(u)
          i=i+1
        }
        
        #melt
        library(reshape)
        Outlier_tmp <- features[,-grep("_sample",colnames(features))]
        Outlier <- Outlier_tmp[,-grep("_ratio",colnames(Outlier_tmp))]
        colnames(Outlier) <- c("Alternative promoter outliers","Expression outliers","Splice outliers",
                               "Chimeric RNA","RNA editing outliers","ASE outliers","SNP",
                               "DNA: Copy-number outlier","DNA: Nucleosome occupancy outliers","DNA: Methylation outliers","UnionN")
        Outlier[is.na(Outlier)] <- 0
        Outlier$Gene <- rownames(Outlier)
        Outlier_plot <- melt(Outlier)
        Outlier_plot$Gene_symbol <- as.character(lapply(strsplit(as.character(Outlier_plot$Gene)," |",fixed = TRUE),function(x) x[1]))
        colnames(Outlier_plot) <- c("Gene","Feature","OutlierN","Gene_symbol")
        
        
        #plot
        #Union
        Union <- Outlier_plot[grep("UnionN",Outlier_plot$Feature),]
        Union$Gene <- factor(Union$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        Union$Gene_symbol <- factor(Union$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        p_Union <- ggplot(Union,aes(x=Gene_symbol,y=OutlierN/sampleN_RNA,fill = Feature))+
          scale_fill_manual(values = c("#666666"))+
          geom_bar(stat = "identity")+
          xlab("")+
          ylab(paste0("RNA+DNA Outliers: ",cutoff,"\nTotal Number: ",sampleN_RNA))+
          guides(fill=FALSE)+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(10,5,-10,5),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            axis.line.y = element_line(color = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            #axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            axis.text.y = element_text(face="bold",  color="black", size=18, angle = 0,hjust=1,vjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0.0","","0.5","","1.0"),expand = c(0,0),limits = c(0,1))
        
        
        #different features contribution
        Features <- Outlier_plot[-grep("UnionN",Outlier_plot$Feature),]
        Features$Gene <- factor(Features$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        Features$Gene_symbol <- factor(Features$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        p_Features <- ggplot(Features,aes(x=Gene_symbol,y=OutlierN/sampleN_RNA,fill=Feature))+
          geom_bar(stat = "identity",position = "fill")+
          xlab("")+
          ylab(paste0("Features: ",cutoff,"\nTotal Number: ",sampleN_RNA))+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(5,5,5,5),units="pt"),
            legend.position="bottom",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            #axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            #axis.text.y = element_text(face="bold",  color="black", size=18, angle = 90,hjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0.01),limits = c(0,1))+
          scale_fill_simpsons(alpha = 0.8)
        
        library(ggpubr)
        ggarrange(p_Union, p_Features,
                  ncol = 1, nrow = 2,heights = c(5,5),align = c("v"))
      }  
      
      #summarize outliers (RNA)
      {
        #input features
        Alt_cutoff <- Alt_targetgenes[grep(cutoff,colnames(Alt_targetgenes))]
        Expression_cutoff <- Expression_targetgenes[grep(cutoff,colnames(Expression_targetgenes))]
        Splicing_cutoff <- Splicing_targetgenes[grep(cutoff,colnames(Splicing_targetgenes))]
        
        #APA_cutoff <- APA_targetgenes[grep(cutoff,colnames(APA_targetgenes))]
        chimeric_cutoff <- chimeric_targetgenes[grep(cutoff,colnames(chimeric_targetgenes))]
        
        Editing_cutoff <- Editing_targetgenes[grep(cutoff,colnames(Editing_targetgenes))]
        ASE_cutoff <- ASE_targetgenes[grep(cutoff,colnames(ASE_targetgenes))]
        SNP_cutoff <- SNP_targetgenes[grep(cutoff,colnames(SNP_targetgenes))]
        
        features <- cbind(Alt_cutoff,Expression_cutoff,Splicing_cutoff,chimeric_cutoff,Editing_cutoff,ASE_cutoff,SNP_cutoff) ## this sequential is important! It determines the label seqeuntial mannually set at melt step!
        features$UnionN <- NA 
        
        #Union number of input features
        i=1
        while(i<=nrow(features)){
          x <- features[i,grep(paste0(cutoff,"_sample"),colnames(features))]
          print(x)
          y <- lapply(x,as.character)
          z <- paste(y,sep = "\n",collapse = "\n")
          u <- unique(unlist(strsplit(z,"\n")))
          u <- u[u!=""]
          u <- u[u!="NA"]
          features[i,which(colnames(features)=="UnionN")] <- length(u)
          i=i+1
        }
        
        #melt
        library(reshape)
        Outlier_tmp <- features[,-grep("_sample",colnames(features))]
        Outlier <- Outlier_tmp[,-grep("_ratio",colnames(Outlier_tmp))]
        colnames(Outlier) <- c("Alternative promoter outliers","Expression outliers","Splice outliers","Chimeric RNA","RNA editing outliers","ASE outliers","SNP","UnionN")
        Outlier[is.na(Outlier)] <- 0
        Outlier$Gene <- rownames(Outlier)
        Outlier_plot <- melt(Outlier)
        Outlier_plot$Gene_symbol <- as.character(lapply(strsplit(as.character(Outlier_plot$Gene)," |",fixed = TRUE),function(x) x[1]))
        colnames(Outlier_plot) <- c("Gene","Feature","OutlierN","Gene_symbol")
        
        
        #plot
        #Union
        RNA_Union <- Outlier_plot[grep("UnionN",Outlier_plot$Feature),]
        RNA_Union$Gene <- factor(Union$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        RNA_Union$Gene_symbol <- factor(Union$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        RNA_Union$background <- NA
        background <- Union
        background$background <- Union$OutlierN
        background$Feature <- "Background"
        background$OutlierN <- NA
        RNA_Union <- rbind(RNA_Union,background)
        p_Union_RNA <- ggplot(RNA_Union)+
          geom_bar(aes(x=Gene_symbol,y=background/sampleN_RNA,fill = Feature),stat = "identity")+
          geom_bar(aes(x=Gene_symbol,y=OutlierN/sampleN_RNA,fill = Feature),stat = "identity")+
          scale_fill_manual(values = c("UnionN"=alpha("black",alpha = 0.5),"Background"=alpha("#666666",alpha = 0.3)))+
          #scale_fill_manual(values = c("UnionN"="#138F6A","Background"=alpha("#666666",alpha = 0.3)))+
          xlab("")+
          ylab(paste0("RNA Outliers: ",cutoff,"\nTotal Number: ",sampleN_RNA))+
          guides(fill=FALSE)+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(10,5,-10,5),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            axis.line.y = element_line(color = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            #axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            axis.text.y = element_text(face="bold",  color="black", size=18, angle = 0,hjust=1,vjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0.0","","0.5","","1.0"),expand = c(0,0),limits = c(0,1))
        
        
        #different features contribution
        RNA_Features <- Outlier_plot[-grep("UnionN",Outlier_plot$Feature),]
        RNA_Features$Gene <- factor(RNA_Features$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        RNA_Features$Gene_symbol <- factor(RNA_Features$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        p_Features_RNA <- ggplot(RNA_Features,aes(x=Gene_symbol,y=OutlierN/sampleN_RNA,fill=Feature))+
          geom_bar(stat = "identity",position = "fill")+
          xlab("")+
          ylab(paste0("Features: ",cutoff,"\nTotal Number: ",sampleN_RNA))+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(5,5,5,5),units="pt"),
            legend.position="bottom",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            #axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            #axis.text.y = element_text(face="bold",  color="black", size=18, angle = 90,hjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0.01),limits = c(0,1))+
          scale_fill_jco(alpha = 0.8)
        
        library(ggpubr)
        ggarrange(p_Union_RNA, p_Features_RNA,
                  ncol = 1, nrow = 2,heights = c(5,5),align = c("v"))
      }
      
      #summarize outliers (DNA)
      {
        #input features
        DNA_CNV_cutoff <- DNA_CNV_targetgenes[grep(cutoff,colnames(DNA_CNV_targetgenes))]
        DNA_CNV_cutoff[,3] <- gsub(".wgs", "", DNA_CNV_cutoff[,3])
        DNA_nucleosome_cutoff <- DNA_nucleosome_targetgenes[grep(cutoff,colnames(DNA_nucleosome_targetgenes))]
        DNA_nucleosome_cutoff[,3] <- gsub(".wgs", "", DNA_nucleosome_cutoff[,3])
        DNA_MeDIP_cutoff <- DNA_MeDIP_targetgenes[grep(cutoff,colnames(DNA_MeDIP_targetgenes))]
        DNA_MeDIP_cutoff[,3] <- gsub(".me", "", DNA_MeDIP_cutoff[,3])
        
        features <- cbind(DNA_CNV_cutoff,DNA_nucleosome_cutoff,DNA_MeDIP_cutoff) ## this sequential is important! It determines the label seqeuntial mannually set at melt step!
        features$UnionN <- NA 
        
        #Union number of input features
        i=1
        while(i<=nrow(features)){
          x <- features[i,grep(paste0(cutoff,"_sample"),colnames(features))]
          y <- lapply(x,as.character)
          z <- paste(y,sep = "\n",collapse = "\n")
          u <- unique(unlist(strsplit(z,"\n")))
          u <- u[u!=""]
          u <- u[u!="NA"]
          features[i,which(colnames(features)=="UnionN")] <- length(u)
          i=i+1
        }
        
        #melt
        library(reshape)
        Outlier_tmp <- features[,-grep("_sample",colnames(features))]
        Outlier <- Outlier_tmp[,-grep("_ratio",colnames(Outlier_tmp))]
        colnames(Outlier) <- c("DNA: Copy-number outlier","DNA: Nucleosome occupancy outliers","DNA: Methylation outliers","UnionN")
        Outlier[is.na(Outlier)] <- 0
        Outlier$Gene <- rownames(Outlier)
        Outlier_plot <- melt(Outlier)
        Outlier_plot$Gene_symbol <- as.character(lapply(strsplit(as.character(Outlier_plot$Gene)," |",fixed = TRUE),function(x) x[1]))
        colnames(Outlier_plot) <- c("Gene","Feature","OutlierN","Gene_symbol")
        
        
        #plot
        #Union
        DNA_Union <- Outlier_plot[grep("UnionN",Outlier_plot$Feature),]
        DNA_Union$Gene <- factor(DNA_Union$Gene,levels=DNA_Union[order(DNA_Union$OutlierN,decreasing = FALSE),]$Gene)
        DNA_Union$Gene_symbol <- factor(DNA_Union$Gene_symbol,levels=DNA_Union[order(DNA_Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        p_Union_DNA <- ggplot(DNA_Union,aes(x=Gene_symbol,y=OutlierN/sampleN_DNA,fill = Feature))+
          scale_fill_manual(values = c("#EEC900"))+
          geom_bar(stat = "identity")+
          xlab("")+
          ylab(paste0("DNA Outliers: ",cutoff,"\nTotal Number: ",sampleN_DNA))+
          guides(fill=FALSE)+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(10,5,-10,5),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.ticks.x = element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            axis.text.x = element_blank(),
            #axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            axis.text.y = element_text(face="bold",  color="black", size=18, angle = 0,hjust=1,vjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0),limits = c(0,1))
        
        #reversed barplot
        {
          DNA_Union$Gene <- factor(DNA_Union$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
          DNA_Union$Gene_symbol <- factor(DNA_Union$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
          DNA_Union$background <- NA
          background <- Union
          background$background <- Union$OutlierN
          background$Feature <- "Background"
          background$OutlierN <- NA
          DNA_Union <- rbind(DNA_Union,background)
          p_Union_DNA_rev <- ggplot(DNA_Union)+
            geom_bar(aes(x=Gene_symbol,y=background/sampleN_DNA,fill = Feature),stat = "identity")+
            geom_bar(aes(x=Gene_symbol,y=OutlierN/sampleN_DNA,fill = Feature),stat = "identity")+
            scale_fill_manual(values = c("UnionN"=alpha("black",alpha = 0.5),"Background"=alpha("#666666",alpha = 0.3)))+
            #scale_fill_manual(values = c("UnionN"="#EEC900","Background"=alpha("#666666",alpha = 0.3)))+
            xlab("")+
            ylab(paste0("DNA Outlier: ",cutoff,"\nTotal Number: ",sampleN_DNA))+
            guides(fill=FALSE)+
            theme_bw()+
            theme(#legend.position="right",
              plot.margin = unit(x=c(5,5,10,5),units="pt"),
              legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
              legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
              plot.title = element_text(hjust = 0.5,size=24,face="bold"),
              axis.ticks.x = element_blank(),
              #axis.text.x = element_blank(),
              axis.line.y = element_line(color = "black"),
              axis.text.x = element_text(color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
              axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
              axis.title.x = element_text(face="bold", color="black", size=20),
              axis.title.y = element_text(face="bold",color="black", size=20))+
            scale_y_reverse(breaks = c(0,0.25,0.50,0.75,1),labels = c("0.0","","0.5","","1.0"),expand = c(0,0),limits = c(1,0))
          }
        
        #different features contribution
        DNA_Features <- Outlier_plot[-grep("UnionN",Outlier_plot$Feature),]
        DNA_Features$Gene <- factor(DNA_Features$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        DNA_Features$Gene_symbol <- factor(DNA_Features$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        p_Features_DNA <- ggplot(DNA_Features,aes(x=Gene_symbol,y=OutlierN/sampleN_DNA,fill=Feature))+
          geom_bar(stat = "identity",position = "fill")+
          xlab("")+
          ylab(paste0("Features: ",cutoff,"\nTotal Number: ",sampleN_DNA))+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(5,5,5,5),units="pt"),
            legend.position="bottom",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            #axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            #axis.text.y = element_text(face="bold",  color="black", size=18, angle = 90,hjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0.01),limits = c(0,1))+
          scale_fill_aaas(alpha = 0.8)
        
        library(ggpubr)
        ggarrange(p_Union_DNA, p_Features_DNA,
                  ncol = 1, nrow = 2,heights = 5,align = c("v"))
      }
      
      #summarize outliers (RNA+DNA)
      {
        #input features
        Alt_cutoff <- Alt_targetgenes[grep(cutoff,colnames(Alt_targetgenes))]
        Alt_cutoff[,3] <- gsub(".pico", "", Alt_cutoff[,3])
        
        Expression_cutoff <- Expression_targetgenes[grep(cutoff,colnames(Expression_targetgenes))]
        Expression_cutoff[,3] <- gsub(".pico", "", Expression_cutoff[,3])
        
        Splicing_cutoff <- Splicing_targetgenes[grep(cutoff,colnames(Splicing_targetgenes))]
        Splicing_cutoff[,3] <- gsub(".pico", "", Splicing_cutoff[,3])
        
        #APA_cutoff <- APA_targetgenes[grep(cutoff,colnames(APA_targetgenes))]
        #APA_cutoff[,3] <- gsub(".pico", "", APA_cutoff[,3])
        
        chimeric_cutoff <- chimeric_targetgenes[grep(cutoff,colnames(chimeric_targetgenes))]
        chimeric_cutoff[,3] <- gsub(".pico", "", chimeric_cutoff[,3])
        
        Editing_cutoff <- Editing_targetgenes[grep(cutoff,colnames(Editing_targetgenes))]
        Editing_cutoff[,3] <- gsub(".pico", "", Editing_cutoff[,3])
        
        ASE_cutoff <- ASE_targetgenes[grep(cutoff,colnames(ASE_targetgenes))]
        ASE_cutoff[,3] <- gsub(".pico", "", ASE_cutoff[,3])
        
        SNP_cutoff <- SNP_targetgenes[grep(cutoff,colnames(SNP_targetgenes))]
        SNP_cutoff[,3] <- gsub(".pico", "", SNP_cutoff[,3])
        
        DNA_CNV_cutoff <- DNA_CNV_targetgenes[grep(cutoff,colnames(DNA_CNV_targetgenes))]
        DNA_CNV_cutoff[,3] <- gsub(".wgs", "", DNA_CNV_cutoff[,3])
        
        DNA_nucleosome_cutoff <- DNA_nucleosome_targetgenes[grep(cutoff,colnames(DNA_nucleosome_targetgenes))]
        DNA_nucleosome_cutoff[,3] <- gsub(".wgs", "", DNA_nucleosome_cutoff[,3])
        
        DNA_MeDIP_cutoff <- DNA_MeDIP_targetgenes[grep(cutoff,colnames(DNA_MeDIP_targetgenes))]
        DNA_MeDIP_cutoff[,3] <- gsub(".me", "", DNA_MeDIP_cutoff[,3])
        
        features <- cbind(DNA_CNV_cutoff,DNA_nucleosome_cutoff,DNA_MeDIP_cutoff,Alt_cutoff,Expression_cutoff,Splicing_cutoff,chimeric_cutoff,Editing_cutoff,ASE_cutoff,SNP_cutoff) ## this sequential is important! It determines the label seqeuntial mannually set at melt step!
        
        #melt
        library(reshape)
        Outlier <- features[,grep("_ratio",colnames(features))]
        colnames(Outlier) <- c("DNA: Copy-number outlier","DNA: Nucleosome occupancy outliers","DNA: Methylation outliers","Alternative promoter outliers","Expression outliers","Splice outliers","Chimeric RNA","RNA editing outliers","ASE outliers","RNA SNP")
        #colnames(Outlier) <- c("DNA_CNV","DNA_nucleosome","DNA_MeDIP","Altpromoter","Expression","Splicing","APA","chimeric","Editing","ASE","SNP")
        Outlier[is.na(Outlier)] <- 0
        Outlier$Gene <- rownames(Outlier)
        Outlier_plot <- melt(Outlier)
        Outlier_plot$Gene_symbol <- as.character(lapply(strsplit(as.character(Outlier_plot$Gene)," |",fixed = TRUE),function(x) x[1]))
        colnames(Outlier_plot) <- c("Gene","Feature","OutlierRatio","Gene_symbol")
        
        
        #plot
        #different features contribution
        Features <- Outlier_plot
        Features$Gene <- factor(Features$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        Features$Gene_symbol <- factor(Features$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        Features_color <- c("DNA: Copy-number outlier"="#FFE600",
                            "DNA: Nucleosome occupancy outliers"="#60AFFE",
                            "DNA: Methylation outliers"="#CDC9C9",
                            "Alternative promoter outliers"="#E9C2A6",
                            "Expression outliers"="#A5435C",
                            "Splice outliers"="#C5E3BF",
                            "Alternative polyadenyltion outliers"="#003F87",
                            "Chimeric RNA"="#FF3D0D",
                            "RNA editing outliers"="#324F17",
                            "ASE outliers"="#800080",
                            "RNA SNP"="#333333")
        p_Features <- ggplot(Features,aes(x=Gene_symbol,y=OutlierRatio,fill=Feature))+
          geom_bar(stat = "identity",position = "fill")+
          xlab("")+
          ylab("")+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(-10,5,-25,5),units="pt"),
            legend.position="bottom",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            #axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            #axis.text.y = element_text(face="bold",  color="black", size=18, angle = 90,hjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0.01),limits = c(0,1))+
          scale_fill_manual(values = Features_color)
        
        library(ggpubr)
        ggarrange(p_Union_RNA,p_Features, p_Union_DNA_rev,
                  ncol = 1, nrow = 3,heights = c(5,2.5,7),align = c("v"), legend = "right", common.legend = TRUE)
      }
    }
    write.csv(Outlier_plot,paste0("./",group,"_Outlier_plot.csv"),quote = FALSE)
    write.csv(DNA_Union,paste0("./",group,"_DNA_Union.csv"),quote = FALSE)
    write.csv(RNA_Union,paste0("./",group,"_RNA_Union.csv"),quote = FALSE)
  }
}

#summary mutation for oncoplot
{
  #Summary mutated genes
  {
    library(data.table)
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/gene_mutation_burden")
    input_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/gene_mutation_burden"
    output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/gene_mutation_burden"
    files <- dir(input_dir)
    files <- grep("VEP",files,value = TRUE)
    
    i=1
    mutation_burden_matrix <- as.data.frame(matrix(numeric(0),ncol=1))
    colnames(mutation_burden_matrix) <- c("Mutated_Gene") 
    mutation_burden_matrix$`Mutated_Gene` <- as.factor(mutation_burden_matrix$`Mutated_Gene`)
    while(i<=length(files)){
      vcf <- read.table(paste0(input_dir,"/",files[i]))
      sample <- as.character(lapply(strsplit(files[i],".",fixed = TRUE),function(x) x[1]))
      colnames(vcf) <- c("Uploaded_variation","Location","Allele","Gene","Feature","Feature_type","Consequence","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","IMPACT","DISTANCE","STRAND","FLAGS","VARIANT_CLASS","SYMBOL","SYMBOL_SOURCE","HGNC_ID","BIOTYPE","CANONICAL","MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM","GENE_PHENO","SIFT","PolyPhen","EXON","INTRON","DOMAINS","miRNA","HGVSc","HGVSp","HGVS_OFFSET","AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF","gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF","MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS")
      #filter COSMIC genes
      vcf_filtered <- vcf[grep("COSV",vcf$Existing_variation),]
      df <- data.frame("Mutated_Gene" = paste(vcf_filtered$Gene,vcf_filtered$SYMBOL,sep="|"),"Variant_Class" = vcf_filtered$VARIANT_CLASS, "Consequence" = vcf_filtered$Consequence)
      mutation_class <- aggregate(Variant_Class ~ Mutated_Gene, df, paste, collapse = ",")
      mutation_class$Variant_Class <- sapply(strsplit(mutation_class$Variant_Class, ","), function(x) paste(rle(x)$values, collapse=";"))
      mutation_class$Variant_Class <- paste0(mutation_class$Variant_Class,";")
      colnames(mutation_class) <- c("Mutated_Gene",sample)
      ##mutation_count <- as.data.frame(setDT(df)[,list(Count=.N),names(df)]) ###
      ##colnames(mutation_count) <- c("Mutated_Gene",sample)
      mutation_burden_matrix <- full_join(mutation_burden_matrix,mutation_class,by=c("Mutated_Gene"="Mutated_Gene"))
      i=i+1
    }
    
    rownames(mutation_burden_matrix) <- mutation_burden_matrix$`Mutated_Gene`
    mutation_burden_matrix <- mutation_burden_matrix[,-which(colnames(mutation_burden_matrix)=="Mutated_Gene")]
    nonsense_row <- which(rownames(mutation_burden_matrix)=="-|-")
    if(length(nonsense_row)!=0){
      mutation_burden_matrix <- mutation_burden_matrix[-which(rownames(mutation_burden_matrix)=="-|-"),]
    } else {
      mutation_burden_matrix <- mutation_burden_matrix
    }
    mutation_burden_matrix <- as.data.frame(t(mutation_burden_matrix))
    
    mutation_burden_matrix$Case.ID <- rownames(mutation_burden_matrix)
    rownames(mutation_burden_matrix) <- seq(from=1,to=nrow(mutation_burden_matrix))
  }
  
  #complexheatmap for mutation
  library(ComplexHeatmap)
  
  
  ##DNA
  {
    #ann_plot <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/oncoplot/RNA_sample_annotation.csv", header = TRUE,)
    #ann_plot$Samples <- gsub("pico","wgs",ann_plot$Samples)
    mat_plot <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/oncoplot/Mutation_burden_class_forplot_DNA.txt",header = TRUE, row.names = 1,sep = "\t")
    mat_plot[] <- lapply(mat_plot, as.character) 
    mat_plot[is.na(mat_plot)] <- ""
    
    mat <- mat_plot[grep("ENSG00000142208|ENSG00000134982|ENSG00000175054|ENSG00000103126|ENSG00000168646|ENSG00000166710|ENSG00000186174|ENSG00000157764|ENSG00000039068|ENSG00000168036|ENSG00000257923|ENSG00000104408|ENSG00000100393|ENSG00000141736|ENSG00000065361|ENSG00000178568|ENSG00000109670|ENSG00000066468|ENSG00000183454|ENSG00000100644|ENSG00000133703|ENSG00000121879|ENSG00000145675|ENSG00000101213|ENSG00000163629|ENSG00000196090|ENSG00000112531|ENSG00000164754|ENSG00000067560|ENSG00000175387|ENSG00000166949|ENSG00000141646|ENSG00000177565|ENSG00000148737|ENSG00000163513|ENSG00000141510|ENSG00000104517|ENSG00000140836",rownames(mat_plot)),]
    #mat <- mat_plot[order,]
    mat <- mat[,-grep("KAPA|NV|Y",colnames(mat))]
    #mat <- mat[,grep("CRC|STAD",colnames(mat))]
    
    #deep sequencing
    mat <- mat[,grep("CRC.PKU.27.wgs|CRC.PKU.28.wgs|CRC.PKU.29.wgs|CRC.PKU.30.wgs|CRC.PKU.32.wgs|NC.PKU.9.wgs|NC.PKU.10.wgs|NC.PKU.11.wgs|NC.PKU.12.wgs|NC.PKU.14.wgs",colnames(mat))]
    #paired sample
    #mat <- mat[,grep("CRC.PKU.3|CRC.PKU.6|CRC.PKU.8|CRC.PKU.9|CRC.PKU.10|CRC.PKU.12|CRC.PKU.13|CRC.PKU.14|CRC.PKU.15|CRC.PKU.16|CRC.PKU.17|CRC.PKU.18|CRC.PKU.21|CRC.PKU.22|CRC.PKU.23|CRC.PKU.24|CRC.PKU.25|CRC.PKU.26|CRC.PKU.27|CRC.PKU.28|CRC.PKU.29|CRC.PKU.30|CRC.PKU.31|CRC.PKU.32|CRC.PKU.33",colnames(mat))]
    
    #mat <- mat[,grep("STAD",colnames(mat))]
    rownames(mat) <- lapply(strsplit(rownames(mat),"|",fixed = TRUE),function(x) x[2])
    mat[is.na(mat)] <- ""
    
    col <- c("SNV" = "#ff7f00", "insertion" = "#984ea3", "deletion" = "#4daf4a", "sequence_alteration" = "#003f87","indel"="#8e2323")
    alter_fun = list(
      background = alter_graphic("rect", fill = "#CCCCCC"),
      SNV = alter_graphic("rect", fill = col["SNV"]),
      insertion = alter_graphic("rect", height = 0.33, fill = col["insertion"]),
      deletion = alter_graphic("rect", height = 0.33, fill = col["deletion"],),
      sequence_alteration = alter_graphic("rect", height = 0.33, fill = col["sequence_alteration"],),
      indel = alter_graphic("rect", height = 0.33, fill = col["indel"],)
    )
    
    column_title <- "Gastritestinal Cancer plasma sample, genes in COSMIC"
    heatmap_legend_param <-
      list(
        title = "Alterations",
        at = c("SNV", "insertion", "deletion","sequence_alteration"),
        labels = c("SNV", "Insertion", "Deletion","Sequence alteration")
      )
    
    oncoPrint(
      mat, alter_fun = alter_fun, col = col,
      #column_order = colnames(mat),
      row_order = unlist(lapply(strsplit(rownames(mat),"|",fixed = TRUE),function(x) x[1])),
      show_column_names = TRUE,
      column_title = column_title,
      heatmap_legend_param = heatmap_legend_param
    )
  }
  ##RNA
  {
    ann_plot <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/oncoplot/RNA_sample_annotation.csv", header = TRUE,)
    mat_plot_RNA <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/oncoplot/Mutation_burden_class_forplot_RNA.txt",header = TRUE, row.names = 1,sep = "\t")
    mat_plot_RNA[] <- lapply(mat_plot_RNA, as.character) 
    mat_plot_RNA[is.na(mat_plot_RNA)] <- ""
    
    mat_RNA <- mat_plot_RNA[grep("ENSG00000142208|ENSG00000134982|ENSG00000175054|ENSG00000103126|ENSG00000168646|ENSG00000166710|ENSG00000186174|ENSG00000157764|ENSG00000039068|ENSG00000168036|ENSG00000257923|ENSG00000104408|ENSG00000100393|ENSG00000141736|ENSG00000065361|ENSG00000178568|ENSG00000109670|ENSG00000066468|ENSG00000183454|ENSG00000100644|ENSG00000133703|ENSG00000121879|ENSG00000145675|ENSG00000101213|ENSG00000163629|ENSG00000196090|ENSG00000112531|ENSG00000164754|ENSG00000067560|ENSG00000175387|ENSG00000166949|ENSG00000141646|ENSG00000177565|ENSG00000148737|ENSG00000163513|ENSG00000141510|ENSG00000104517|ENSG00000140836",rownames(mat_plot_RNA)),]
    #mat_RNA <- mat_RNA[,grep("STAD|CRC",colnames(mat_RNA))]
    rownames(mat_RNA) <- lapply(strsplit(rownames(mat_RNA),"|",fixed = TRUE),function(x) x[2])
    mat_RNA[is.na(mat_RNA)] <- ""
    
    #deep sequencing
    mat_RNA <- mat_RNA[,grep("CRC.PKU.27.pico|CRC.PKU.28.pico|CRC.PKU.29.pico|CRC.PKU.30.pico|CRC.PKU.32.pico|NC.PKU.9.pico|NC.PKU.10.pico|NC.PKU.11.pico|NC.PKU.12.pico|NC.PKU.13.pico|NC.PKU.14.pico",colnames(mat_RNA))]
    #mat_RNA <- mat_RNA[,grep("STAD|CRC",colnames(mat_RNA))]
    #paired sample
    
    #ann <- ann_plot
    ann <- ann_plot[grep("STAD|CRC",ann_plot$Samples),]
    
    #deep sequencing
    #ann <- ann[grep("CRC-PKU-27-pico|CRC-PKU-28-pico|CRC-PKU-29-pico|CRC-PKU-30-pico|CRC-PKU-32-pico|NC-PKU-9-pico|NC-PKU-10-pico|NC-PKU-11-pico|NC-PKU-12-pico|NC-PKU-14-pico",ann$Samples),]
    
    ann <- ann[,-which(colnames(ann)=="Samples")]
    
    
    
    
    colours_ann <- list(
      'Gender'=c('Male'='#5CACEE','Female'='#FF0000'),
      'Stage'=c('Stage I'='#DCDCDC','Stage II'='#A3A3A3','Stage III'='#4A4A4A','Stage IV'='#000000'),
      'Location'=c('No biopsy'='#FFFFFF','Others'='#DCDCDC',
                   'Rectum'='#FFF8DC','Sigmoid colon'='#EEDD82','Descending colon'='#FEE5AC','Transverse colon'='#EDCB62','Ascending colon'='#FCB514',
                   'Fundus'='#FEF1B5','Body'='#E5BC3B','Antrum'='#CD9B1D'),
      'MMR'=c('No biopsy'='#FFFFFF','p'='#EDCB62','d'='#FCB514'),
      'HER2'=c('No biopsy'='#FFFFFF','-'='#E5BC3B','+'='#CD9B1D'),
      'Type'=c('Colorectal cancer'='#FCB514','Stomach cancer'='#CD9B1D')
    )
    colAnn <- HeatmapAnnotation(col = colours_ann,
                                #cbar = anno_oncoprint_barplot(),
                                df = ann,which = 'col')
    
    
    col <- c("SNV" = "#ff7f00", "insertion" = "#984ea3", "deletion" = "#4daf4a", "sequence_alteration" = "#003f87","indel"="#8e2323")
    alter_fun = list(
      background = alter_graphic("rect", fill = "#CCCCCC"),
      SNV = alter_graphic("rect", fill = col["SNV"]),
      insertion = alter_graphic("rect", height = 0.33, fill = col["insertion"]),
      deletion = alter_graphic("rect", height = 0.33, fill = col["deletion"],),
      sequence_alteration = alter_graphic("rect", height = 0.33, fill = col["sequence_alteration"],),
      indel = alter_graphic("rect", height = 0.33, fill = col["indel"],)
    )
    
    column_title <- "Gastrointestinal cancer plasma sample, genes in COSMIC"
    heatmap_legend_param <-
      list(
        title = "Alterations",
        at = c("SNV", "insertion", "deletion","sequence_alteration","indel"),
        labels = c("SNV", "Insertion", "Deletion","Sequence alteration","Indel")
      )
    
    #simple
    oncoPrint(
      mat_RNA, alter_fun = alter_fun, col = col,
      #column_order = colnames(mat),
      row_order = unlist(lapply(strsplit(rownames(mat_RNA),"|",fixed = TRUE),function(x) x[1])),
      show_column_names = TRUE,
      column_title = column_title,
      heatmap_legend_param = heatmap_legend_param
    )
    
    oncoPrint(
      mat_RNA, alter_fun = alter_fun, col = col,
      #column_order = colnames(mat),
      #row_order = unlist(lapply(strsplit(order,"|",fixed = TRUE),function(x) x[2])),
      show_column_names = FALSE,
      column_title = column_title,
      heatmap_legend_param = heatmap_legend_param,
      bottom_annotation = colAnn,
      top_annotation = HeatmapAnnotation(
        cbar = anno_oncoprint_barplot(),
        show_annotation_name = TRUE)
      #  Age = 1:ncol(mat_RNA),
      #  Gender = 1:ncol(mat_RNA),
      #  Stage = 1:ncol(mat_RNA)
      #)
    )
  }
}

#concordant mutation
{
#demo
{
  VEP <- read.table("./VEP_RNA/CRC-PKU-27-pico.rmEDIT.filtered.VEP.txt",seq = "\t")
  concordent <- read.table("../individual_compare/CRC-PKU-27.concordant.vcf.diff.sites_in_files",sep = "\t",header = TRUE)
  concordent <- concordent[which(concordent$IN_FILE=="B"),]
  SNPs <- paste0(concordent$CHROM,"_",concordent$POS1,"_",concordent$REF1,"/",concordent$ALT1)
  SNPs <- paste0(concordent$CHROM,":",concordent$POS1)
    
  i=1
  concordent_VEP={}
  while(i<=length(SNPs)){
  SNPs_VEP <- VEP[grep(SNPs[i],VEP$V2),]
  concordent_VEP <- rbind(concordent_VEP,SNPs_VEP)
  i=i+1
  }
}

#oncoplot
{
 concordent_mutation <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/oncoplot/Mutation_burden_class_forplot_concordent_Mutation.txt",sep = "\t",header = TRUE)
  CRC <- apply(concordent_mutation[,grep("CRC",colnames(concordent_mutation))], 1,function(x) length(na.omit(x)))
  CRC_filtered <- CRC[CRC>2]
  CRC_concordent_mutation_filtered <- concordent_mutation[names(CRC_filtered),grep("CRC",colnames(concordent_mutation))]
  CRC_concordent_mutation_filtered[] <- lapply(CRC_concordent_mutation_filtered, as.character) 
  CRC_concordent_mutation_filtered[is.na(CRC_concordent_mutation_filtered)] <- ""
  rownames(CRC_concordent_mutation_filtered) <- as.character(lapply(strsplit(rownames(CRC_concordent_mutation_filtered),"\\|"),function(x) x[2]))
  colnames(CRC_concordent_mutation_filtered) <- gsub(".","-",fixed = TRUE,colnames(CRC_concordent_mutation_filtered))
  
  NC <- apply(concordent_mutation[,grep("NC",colnames(concordent_mutation))], 1,function(x) length(na.omit(x)))
  NC_filtered <- NC[NC>2]
  NC_concordent_mutation_filtered <- concordent_mutation[names(NC_filtered),grep("NC",colnames(concordent_mutation))]
  NC_concordent_mutation_filtered[] <- lapply(NC_concordent_mutation_filtered, as.character) 
  NC_concordent_mutation_filtered[is.na(NC_concordent_mutation_filtered)] <- ""
  rownames(NC_concordent_mutation_filtered) <- as.character(lapply(strsplit(rownames(NC_concordent_mutation_filtered),"\\|"),function(x) x[2]))
  NC_concordent_mutation_filtered <- NC_concordent_mutation_filtered[,c(6,1,2,3,4,5)]
  colnames(NC_concordent_mutation_filtered) <- gsub(".","-",fixed = TRUE,colnames(NC_concordent_mutation_filtered))
  
  col <- c("SNV" = "#ff7f00", "insertion" = "#984ea3", "deletion" = "#4daf4a", "sequence_alteration" = "#003f87","indel"="#af4035")
  alter_fun = list(
    background = alter_graphic("rect", fill = "white"),#CCCCCC
    SNV = alter_graphic("rect", fill = col["SNV"]),
    insertion = alter_graphic("rect", height = 0.33, fill = col["insertion"]),
    deletion = alter_graphic("rect", height = 0.33, fill = col["deletion"],),
    sequence_alteration = alter_graphic("rect", height = 0.33, fill = col["sequence_alteration"],),
    indel = alter_graphic("rect", height = 0.33, fill = col["indel"],)
  )
  
  column_title <- ""
  heatmap_legend_param <-
    list(
      title = "Alterations",
      at = c("SNV", "insertion", "deletion","sequence_alteration","indel"),
      labels = c("SNV", "Insertion", "Deletion","Sequence alteration","Indel")
    )
  
  #NC
  oncoPrint(
    NC_concordent_mutation_filtered, alter_fun = alter_fun, col = col,
    left_annotation = NULL,
    right_annotation = NULL,
    bottom_annotation = NULL,
    top_annotation = NULL,
    column_order = colnames(NC_concordent_mutation_filtered),
    remove_empty_columns = TRUE,
    #row_order = unlist(lapply(strsplit(order,"|",fixed = TRUE),function(x) x[2])),
    show_column_names = TRUE,
    column_title = column_title,
    heatmap_legend_param = heatmap_legend_param
  )
  
  #CRC
  oncoPrint(
    CRC_concordent_mutation_filtered, alter_fun = alter_fun, col = col,
    left_annotation = NULL,
    right_annotation = NULL,
    bottom_annotation = NULL,
    top_annotation = NULL,
    column_order = colnames(CRC_concordent_mutation_filtered),
    remove_empty_columns = TRUE,
    #row_order = unlist(lapply(strsplit(order,"|",fixed = TRUE),function(x) x[2])),
    show_column_names = TRUE,
    column_title = column_title,
    heatmap_legend_param = heatmap_legend_param
  )
}
  #function of concordant mutated genes
  genes <- names(NC[NC>1])
  output_KEGG <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/Function of concordently mutated gene/Mutation_NC_KEGG.txt"
  output_GO_BP <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/Function of concordently mutated gene/Mutation_NC_GO_BP.txt"
  output_GO_CC <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/Function of concordently mutated gene/Mutation_NC_GO_CC.txt"
  output_GO_MF <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/Function of concordently mutated gene/Mutation_NC_GO_MF.txt"
  KEGG_GO(genes,output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC)
  
  genes <- names(CRC[CRC>1])
  output_KEGG <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/Function of concordently mutated gene/Mutation_CRC_KEGG.txt"
  output_GO_BP <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/Function of concordently mutated gene/Mutation_CRC_GO_BP.txt"
  output_GO_CC <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/Function of concordently mutated gene/Mutation_CRC_GO_CC.txt"
  output_GO_MF <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/Function of concordently mutated gene/Mutation_CRC_GO_MF.txt"
  KEGG_GO(genes,output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC)
}