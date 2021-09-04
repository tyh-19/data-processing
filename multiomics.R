##multiomics

#library preparation
{
library(edgeR)
library(ggplot2) #plot
library(sunburstR) #sunburst plot
library(htmlwidgets) #for html sunburst plot
library(extrafont) #for plot font load
fonts() #load fonts
library(ggsci) #for plot color theme 
library(cowplot) #for plot arrange
library(tidyverse) #for dataframe process
library(progress) #for progress bar
}

#published data number
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/01.Sample summary")
  
  #RNA
  published <- read.csv("published_RNA.csv",header = TRUE)
  published <- published[-grep("Controls-redundant",published$Disease_Type),]
  published <- published[-grep("Controls-non-cancer",published$Disease_Type),]
  published <- published[-grep("Prostate cancer",published$Disease_Type),]
  published <- published[-grep("Renal cancer",published$Disease_Type),]
  published <- published[-grep("Cholangio carcinoma",published$Disease_Type),]
  published <- published[-grep("Pregnant",published$Disease_Type),]
  published <- published[-grep("Spontaneous preterm birth",published$Disease_Type),]
  published <- published[-grep("Alzheimers disease",published$Disease_Type),]
  published <- published[-grep("Coronary heart disease",published$Disease_Type),]
  published <- published[-grep("Urine",published$Specimen),]
  published$Disease_Type <- factor(published$Disease_Type,levels = unique(published[order(published$Sum_size,decreasing= TRUE),]$`Disease_Type`))
  #published$Disease_Type <- factor(published$Disease_Type,levels = c("Controls","Lung cancer","Pancreatic cancer","Liver cancer","Breast cancer","Colorectal cancer","Glioblastoma","Stomach cancer","Esophageal cancer","Cholangio carcinoma"))
  published$RNA_Source <- factor(published$RNA_Source,levels = c("Cellular","Cell-free","EV","Platelet"))
  p_RNA <- ggplot(published,aes(x=Disease_Type,y=Sample_Size,fill = RNA_Source))+geom_bar(stat = "identity")+
    geom_text(aes(x=Disease_Type, y=Sum_size+50,label=Sum_size),size = 6,angle = 90)+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(5,5,-10,5),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_blank(),
      #axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 0,vjust = 0.5),
      axis.text.y = element_text(face="bold",  color="black", size=18, angle = 90,hjust=0.5),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    xlab("")+
    scale_fill_jco(alpha = 0.8)
  
  #DNA
  published_DNA <- read.csv("published_DNA.csv",header = TRUE)
  published_DNA <- published_DNA[-grep("Liver disease",published_DNA$Disease_Type),]
  published_DNA <- published_DNA[-grep("Leukemia",published_DNA$Disease_Type),]
  published_DNA <- published_DNA[-grep("Pediatric medullablastoma",published_DNA$Disease_Type),]
  published_DNA <- published_DNA[-grep("Prostate cancer",published_DNA$Disease_Type),]
  published_DNA <- published_DNA[-grep("Colonrectal disease",published_DNA$Disease_Type),]
  published_DNA <- published_DNA[-grep("Pancreas disease",published_DNA$Disease_Type),]
  published_DNA <- published_DNA[-grep("Stomach disease",published_DNA$Disease_Type),]
  published_DNA <- published_DNA[-grep("Thyroid disease",published_DNA$Disease_Type),]
  published_DNA <- published_DNA[-grep("Thyroid cancer",published_DNA$Disease_Type),]
  published_DNA$Disease_Type <- factor(published_DNA$Disease_Type,levels = levels(published$Disease_Type))
  #published_DNA$Disease_Type <- factor(published_DNA$Disease_Type,levels = unique(published_DNA[order(published_DNA$Sum_size,decreasing= TRUE),]$`Disease_Type`))
  published_DNA$Library_Type <- factor(published_DNA$Library_Type,levels = c("WGS","target DNA-seq","5mC WGBS","5mC DIP-seq","5hmC DIP-seq","5mC BS-seq"))
  p_DNA <- ggplot(published_DNA,aes(x=Disease_Type,y=Sample_Size,fill = Library_Type))+geom_bar(stat = "identity")+
    geom_text(aes(x=Disease_Type, y=Sum_size+180,label=Sum_size),size = 6,angle = 90)+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-20,5,10,5),units="pt"),
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      legend.position="right",
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(face="bold",  color="black", size=18,angle = 90,hjust=0.5),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    scale_y_reverse()+xlab("")+scale_fill_aaas(alpha = 0.7)+
    scale_x_discrete(position = "top")
  
  ggarrange(p_RNA, p_DNA,
            ncol = 1, nrow = 2,heights = 5,align = c("v"))
}

#published data: differential expression between CRC and NC (Tumor and Normal) in plasma/PBMC/tissue, rank by foldchange genes (FDR<0.05) and then enrich at pathway level by GSEA
{
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/09.Published_data/")
    #GSE142978_HCC_vs_NC
    mat_raw <- read.table("./expression_matrix/GSE142987_sample_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    des <- read.csv("./group/GSE142978_HCC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE142978_HCC_vs_NC"
    #GSE174302
    mat_raw <- read.table("./expression_matrix/GSE174302_intron-spanning-available-samples.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    ##CRCvsNC
    des <- read.csv("./group/GSE174302_CRC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE174302_CRC_vs_NC"
    ##ESCAvsNC
    des <- read.csv("./group/GSE174302_ESCA_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE174302_ESCA_vs_NC"
    ##HCCvsNC
    des <- read.csv("./group/GSE174302_HCC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE174302_HCC_vs_NC"
    ##LUADvsNC
    des <- read.csv("./group/GSE174302_LUAD_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE174302_LUAD_vs_NC"
    ##STADvsNC
    des <- read.csv("./group/GSE174302_STAD_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE174302_STAD_vs_NC"
    #GSE68086
    mat_raw <- read.table("./expression_matrix/GSE68086_TEP_data_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    ##BRCAvsNC
    des <- read.csv("./group/GSE68086_BRCA_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE68086_BRCA_vs_NC"
    ##CRCvsNC
    des <- read.csv("./group/GSE68086_CRC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE68086_CRC_vs_NC"
    ##GBMvsNC
    des <- read.csv("./group/GSE68086_GBM_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE68086_GBM_vs_NC"
    ##HepatobiliaryvsNC
    des <- read.csv("./group/GSE68086_Hepatobiliary_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE68086_Hepatobiliary_vs_NC"
    ##NSCLCvsNC
    des <- read.csv("./group/GSE68086_NSCLC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE68086_NSCLC_vs_NC"
    ##PDACvsNC
    des <- read.csv("./group/GSE68086_PDAC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE68086_PDAC_vs_NC"
    #GSE89843
    mat_raw <- read.table("./expression_matrix/GSE89843_TEP_Count_Matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    ##Multipl-SclerosisvsNC
    des <- read.csv("./group/GSE89843_Multipl-Sclerosis_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE89843_Multiple-Sclerosis_vs_NC"
    ##Pulmonary-HypertensionvsNC
    des <- read.csv("./group/GSE89843_Pulmonary-Hypertension_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE89843_Pulmonary-Hypertension_vs_NC"
    ##Chronic-PancreatitisvsNC
    des <- read.csv("./group/GSE89843_Chronic-Pancreatitis_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE89843_Chronic-Pancreatitis_vs_NC"
    ##EpilepsyvsNC
    des <- read.csv("./group/GSE89843_Epilepsy_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE89843_Epilepsy_vs_NC"
    ##Non-significant-AtherosclerosisvsNC
    des <- read.csv("./group/GSE89843_Atherosclerosis_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE89843_Atherosclerosis_vs_NC"
    ##Angina-PectorisvsNC
    des <- read.csv("./group/GSE89843_Angina-Pectoris_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE89843_Angina-Pectoris_vs_NC"
    ##NSCLCvsNC
    des <- read.csv("./group/GSE89843_NSCLC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE89843_NSCLC_vs_NC"
    #exoRBase
    mat_raw <- read.table("./expression_matrix/exoRbase_featurecounts_inhouse.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    ##HCCvsNC
    des <- read.csv("./group/exoRBase_HCC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "exoRBase_HCC_vs_NC"
    ##CRCvsNC
    des <- read.csv("./group/exoRBase_CRC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "exoRBase_CRC_vs_NC"
    ##CHDvsNC
    des <- read.csv("./group/exoRBase_CHD_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "exoRBase_CHD_vs_NC"
    ##PAADvsNC
    des <- read.csv("./group/exoRBase_PAAD_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "exoRBase_PAAD_vs_NC"
    #GSE133684
    mat_raw_forward <- read.table("./expression_matrix/GSE133684_featurecounts_forward_inhouse.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    mat_raw_reverse <- read.table("./expression_matrix/GSE133684_featurecounts_reverse_inhouse.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    mat_raw_forward$gene_id <- rownames(mat_raw_forward)
    mat_raw_reverse$gene_id <- rownames(mat_raw_reverse)
    mat_raw <- join(mat_raw_forward,mat_raw_reverse, by = "gene_id", type = "full")
    rownames(mat_raw) <- mat_raw$gene_id
    mat_raw[is.na(mat_raw)] <- 0
    mat_raw <- mat_raw[,-grep("gene_id",colnames(mat_raw))]
    ##PDACvsNC
    des <- read.csv("./group/GSE133684_PDAC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE133684_PDAC_vs_NC"
    #GSE40174
    mat_raw <- read.table("./expression_matrix/GSE40174_featurecounts_inhouse.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    ##PDACvsNC
    des <- read.csv("./group/GSE40174_PDAC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE40174_PDAC_vs_NC"
    #GSE120663
    mat_raw <- read.table("./expression_matrix/GSE120663_featurecounts_inhouse.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    ##HCCvsNC
    des <- read.csv("./group/GSE120663_HCC_vs_NC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE120663_HCC_vs_NC"
    
    #GSE117623
    mat_raw <- read.table("./expression_matrix/GSE117623_featurecounts_inhouse.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    ##HCCvsNC
    des <- read.csv("./group/GSE117623_HCC_vs_DC.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "GSE117623_HCC_vs_DC"
    
    #TCGA_COAD
    mat_raw <- read.table("./expression_matrix/TCGA_COAD_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    des <- read.csv("./group/TCGA_COAD_PrimaryTumor_vs_TissueNormal.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "TCGA_COAD_PrimaryTumor_vs_TissueNormal"
    
    #TCGA_STAD
    mat_raw <- read.table("./expression_matrix/TCGA_STAD_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    des <- read.csv("./group/TCGA_STAD_PrimaryTumor_vs_TissueNormal.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "TCGA_STAD_PrimaryTumor_vs_TissueNormal"
    
    #TCGA_READ
    mat_raw <- read.table("./expression_matrix/TCGA_READ_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    des <- read.csv("./group/TCGA_READ_PrimaryTumor_vs_TissueNormal.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "TCGA_READ_PrimaryTumor_vs_TissueNormal"
    
    #TCGA_BRCA
    mat_raw <- read.table("./expression_matrix/TCGA_BRCA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    ##all
    des <- read.csv("./group/TCGA_BRCA_PrimaryTumor_vs_TissueNormal.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "TCGA_BRCA_PrimaryTumor_vs_TissueNormal"
    ##paired
    des <- read.csv("./group/TCGA_BRCA_paired_PrimaryTumor_vs_TissueNormal.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "TCGA_BRCA_paired_PrimaryTumor_vs_TissueNormal"
    
    #TCGA_ESCA
    mat_raw <- read.table("./expression_matrix/TCGA_ESCA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    des <- read.csv("./group/TCGA_ESCA_PrimaryTumor_vs_TissueNormal.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "TCGA_ESCA_PrimaryTumor_vs_TissueNormal"
    
    #TCGA_LIHC
    mat_raw <- read.table("./expression_matrix/TCGA_LIHC_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    des <- read.csv("./group/TCGA_LIHC_PrimaryTumor_vs_TissueNormal.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "TCGA_LIHC_PrimaryTumor_vs_TissueNormal"
    
    #TCGA_LUAD
    mat_raw <- read.table("./expression_matrix/TCGA_LUAD_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    des <- read.csv("./group/TCGA_LUAD_PrimaryTumor_vs_TissueNormal.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "TCGA_LUAD_PrimaryTumor_vs_TissueNormal"
    
    #TCGA_LUSC
    mat_raw <- read.table("./expression_matrix/TCGA_LUSC_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    des <- read.csv("./group/TCGA_LUSC_PrimaryTumor_vs_TissueNormal.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "TCGA_LUSC_PrimaryTumor_vs_TissueNormal"
    
    #TCGA_THCA
    mat_raw <- read.table("./expression_matrix/TCGA_THCA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    des <- read.csv("./group/TCGA_THCA_PrimaryTumor_vs_TissueNormal.csv", header = TRUE, check.names=FALSE, sep=',')
    dataset <- "TCGA_THCA_PrimaryTumor_vs_TissueNormal"
  }
  {  
    samples <- des$samples
    #samples <- gsub(".","-",samples,fixed=TRUE)
    #samples <- gsub("X","",samples,fixed=TRUE)
    group <- des$group
    #batch <- des$batch
    
    i=1
    mat=as.data.frame(array(dim=c(length(rownames(mat_raw)),1)))
    while (i<=length(samples)) {
      temp <- mat_raw[,which(colnames(mat_raw)==samples[i])]
      #temp <- mat_raw[,grep(samples[i],colnames(mat_raw),fixed=TRUE)]
      #print(i)
      #print(which(colnames(mat_raw)==samples[i]))
      temp <- as.data.frame(temp)
      colnames(temp) <- samples[i]
      mat <- cbind(mat,temp)
      i=i+1
    }
    
    mat <- mat[,-1]
    rownames(mat) <- rownames(mat_raw)
    
    y <- DGEList(counts=mat, samples=samples, group=group)     #规定输入格式，必须包含这几个参数
    y <- calcNormFactors(y, method="TMM")                      #归一化处理
    y <- estimateDisp(y)
    test <- exactTest(y, pair = c("negative","positive"), dispersion = "auto")
    #https://support.bioconductor.org/p/64807/
    #广义线性模型可以针对多因素的情况，在这个场景下，只比较cancer与normal用exactTest就足够了
    {
      #design <- model.matrix(~0+group)   
      #y <- estimateDisp(y, design)          #广义线性模型计算离散度（common&trended&tagwise）
      #fit <- glmFit(y, design)              #差异表达分析函数                     
      #test <- glmLRT(fit, coef=2)           #差异表达分析函数 
    }
    
    res <- topTags(test, n=nrow(mat), sort.by='none')     #输出计算差异基因的结果
    res <- cbind(res$table, baseMean=2^(res$table$logCPM)) #输出差异基因结果加上baseMean一栏
    mapped_names <- colnames(res)                   #后面都是整理+输出
    for(i in 1:ncol(res)){
      if(colnames(res)[i] == 'logFC'){
        mapped_names[i] <- 'log2FoldChange'
      }else if(colnames(res)[i] == 'PValue'){
        mapped_names[i] <- 'pvalue'
      }else if(colnames(res)[i] == 'FDR') {
        mapped_names[i] <- 'padj'
      }else{
        mapped_names[i] <- colnames(res)[i]
      }
    }
    colnames(res) <- mapped_names
    write.table(res, paste0("./output/",dataset,"_edger_exact.csv"), sep=',', quote=FALSE, row.names=TRUE) #输出文件，规定名字
  }
  
  #make ranked genelist
  {
    res.sort <- res[res$log2FoldChange!=0,]
    res.sort <- res.sort[sort(res.sort$log2FoldChange,decreasing = T, index.return = T)$ix,]
    res.sort <- res.sort[grep("ENSG",rownames(res.sort)),]
    
    genelist <- as.numeric(res.sort$log2FoldChange)
    names(genelist) <- as.character(lapply(strsplit(rownames(res.sort),"\\|"), function(x) x[1]))
    names(genelist) <- as.character(lapply(strsplit(names(genelist),".",fixed = TRUE), function(x) x[1]))
    
  }
  
  #make geneset from KEGG
  {
    PATH_ID_NAME <- read.csv("./pathway/PATH_ID_NAME_KEGGplusHallmark.txt",header = TRUE, row.names = 1,sep = "\t",check.names=FALSE)
    #gmt <- PATH_ID_NAME[,c(4,1)] # local_modified
    gmt <- PATH_ID_NAME[,c(6,1)] # local_plusHallmarker
    colnames(gmt) <- c("ont","gene")
  }
  
  #GSEA by clusterprofile and plot
  {
    egmt2 <- GSEA(genelist,minGSSize = 0, maxGSSize = 1000,
                  pvalueCutoff = 1, TERM2GENE = gmt)
    write.table(egmt2@result,paste0("./output/",dataset,"_GSEA.txt"),quote = FALSE,sep = "\t")
    
    #gsea <- gseaplot2(egmt2,1)
    #pdf("Ribosome_PBMC_CRCvsNC.pdf")
    #plot(gsea)
    #dev.off()
  }
  
  #bubble plot
  {
    GSE142978_HCC_vs_NC <- read.csv("./output/GSE142978_HCC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE142978_HCC_vs_NC$dataset <- "GSE142978_HCC_vs_NC"

    GSE174302_CRC_vs_NC <- read.csv("./output/GSE174302_CRC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE174302_CRC_vs_NC$dataset <- "GSE174302_CRC_vs_NC"
    
    GSE174302_ESCA_vs_NC <- read.csv("./output/GSE174302_ESCA_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE174302_ESCA_vs_NC$dataset <- "GSE174302_ESCA_vs_NC"
    
    GSE174302_HCC_vs_NC <- read.csv("./output/GSE174302_HCC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE174302_HCC_vs_NC$dataset <- "GSE174302_HCC_vs_NC"
    
    GSE174302_LUAD_vs_NC <- read.csv("./output/GSE174302_LUAD_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE174302_LUAD_vs_NC$dataset <- "GSE174302_LUAD_vs_NC"
    
    GSE174302_STAD_vs_NC <- read.csv("./output/GSE174302_STAD_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE174302_STAD_vs_NC$dataset <- "GSE174302_STAD_vs_NC"
    
    GSE68086_BRCA_vs_NC <- read.csv("./output/GSE68086_BRCA_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE68086_BRCA_vs_NC$dataset <- "GSE68086_BRCA_vs_NC"
    
    GSE68086_CRC_vs_NC <- read.csv("./output/GSE68086_CRC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE68086_CRC_vs_NC$dataset <- "GSE68086_CRC_vs_NC"
    
    GSE68086_GBM_vs_NC <- read.csv("./output/GSE68086_GBM_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE68086_GBM_vs_NC$dataset <- "GSE68086_GBM_vs_NC"
    
    GSE68086_Hepatobiliary_vs_NC <- read.csv("./output/GSE68086_Hepatobiliary_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE68086_Hepatobiliary_vs_NC$dataset <- "GSE68086_Hepatobiliary_vs_NC"
    
    GSE68086_NSCLC_vs_NC <- read.csv("./output/GSE68086_NSCLC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE68086_NSCLC_vs_NC$dataset <- "GSE68086_NSCLC_vs_NC"
    
    GSE68086_PDAC_vs_NC <- read.csv("./output/GSE68086_PDAC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE68086_PDAC_vs_NC$dataset <- "GSE68086_PDAC_vs_NC"
    
    GSE89843_Multiple_Sclerosis_vs_NC <- read.csv("./output/GSE89843_Multiple-Sclerosis_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE89843_Multiple_Sclerosis_vs_NC$dataset <- "GSE89843_Multiple-Sclerosis_vs_NC"
    
    GSE89843_Pulmonary_Hypertension_vs_NC <- read.csv("./output/GSE89843_Pulmonary-Hypertension_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE89843_Pulmonary_Hypertension_vs_NC$dataset <- "GSE89843_Pulmonary-Hypertension_vs_NC"
    
    GSE89843_Chronic_Pancreatitis_vs_NC <- read.csv("./output/GSE89843_Chronic-Pancreatitis_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE89843_Chronic_Pancreatitis_vs_NC$dataset <- "GSE89843_Chronic-Pancreatitis_vs_NC"
    
    GSE89843_Epilepsy_vs_NC <- read.csv("./output/GSE89843_Epilepsy_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE89843_Epilepsy_vs_NC$dataset <- "GSE89843_Epilepsy_vs_NC"
    
    GSE89843_Atherosclerosis_vs_NC <- read.csv("./output/GSE89843_Atherosclerosis_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE89843_Atherosclerosis_vs_NC$dataset <- "GSE89843_Atherosclerosis_vs_NC"
    
    GSE89843_Angina_Pectoris_vs_NC <- read.csv("./output/GSE89843_Angina-Pectoris_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE89843_Angina_Pectoris_vs_NC$dataset <- "GSE89843_Angina-Pectoris_vs_NC"
    
    GSE89843_NSCLC_vs_NC <- read.csv("./output/GSE89843_NSCLC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE89843_NSCLC_vs_NC$dataset <- "GSE89843_NSCLC_vs_NC"
    
    exoRBase_HCC_vs_NC <- read.csv("./output/exoRBase_HCC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    exoRBase_HCC_vs_NC$dataset <- "exoRBase_HCC_vs_NC"
    
    exoRBase_CRC_vs_NC <- read.csv("./output/exoRBase_CRC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    exoRBase_CRC_vs_NC$dataset <- "exoRBase_CRC_vs_NC"
    
    exoRBase_CHD_vs_NC <- read.csv("./output/exoRBase_CHD_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    exoRBase_CHD_vs_NC$dataset <- "exoRBase_CHD_vs_NC"
    
    exoRBase_PAAD_vs_NC <- read.csv("./output/exoRBase_PAAD_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    exoRBase_PAAD_vs_NC$dataset <- "exoRBase_PAAD_vs_NC"
    
    GSE133684_PDAC_vs_NC <- read.csv("./output/GSE133684_PDAC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE133684_PDAC_vs_NC$dataset <- "GSE133684_PDAC_vs_NC"
    
    GSE40174_PDAC_vs_NC <- read.csv("./output/GSE40174_PDAC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE40174_PDAC_vs_NC$dataset <- "GSE40174_PDAC_vs_NC"
    
    GSE120663_HCC_vs_NC <- read.csv("./output/GSE120663_HCC_vs_NC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE120663_HCC_vs_NC$dataset <- "GSE120663_HCC_vs_NC"
    
    GSE117623_HCC_vs_DC <- read.csv("./output/GSE117623_HCC_vs_DC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    GSE117623_HCC_vs_DC$dataset <- "GSE117623_HCC_vs_DC"
    
    TCGA_COAD_PrimaryTumor_vs_TissueNormal <- read.csv("./output/TCGA_COAD_PrimaryTumor_vs_TissueNormal_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    TCGA_COAD_PrimaryTumor_vs_TissueNormal$dataset <- "TCGA_COAD_PrimaryTumor_vs_TissueNormal"
    
    TCGA_STAD_PrimaryTumor_vs_TissueNormal <- read.csv("./output/TCGA_STAD_PrimaryTumor_vs_TissueNormal_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    TCGA_STAD_PrimaryTumor_vs_TissueNormal$dataset <- "TCGA_STAD_PrimaryTumor_vs_TissueNormal"
    
    TCGA_BRCA_PrimaryTumor_vs_TissueNormal <- read.csv("./output/TCGA_BRCA_PrimaryTumor_vs_TissueNormal_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    TCGA_BRCA_PrimaryTumor_vs_TissueNormal$dataset <- "TCGA_BRCA_PrimaryTumor_vs_TissueNormal"
    
    TCGA_ESCA_PrimaryTumor_vs_TissueNormal <- read.csv("./output/TCGA_ESCA_PrimaryTumor_vs_TissueNormal_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    TCGA_ESCA_PrimaryTumor_vs_TissueNormal$dataset <- "TCGA_ESCA_PrimaryTumor_vs_TissueNormal"
    
    TCGA_LIHC_PrimaryTumor_vs_TissueNormal <- read.csv("./output/TCGA_LIHC_PrimaryTumor_vs_TissueNormal_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    TCGA_LIHC_PrimaryTumor_vs_TissueNormal$dataset <- "TCGA_LIHC_PrimaryTumor_vs_TissueNormal"
    
    TCGA_LUAD_PrimaryTumor_vs_TissueNormal <- read.csv("./output/TCGA_LUAD_PrimaryTumor_vs_TissueNormal_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    TCGA_LUAD_PrimaryTumor_vs_TissueNormal$dataset <- "TCGA_LUAD_PrimaryTumor_vs_TissueNormal"
    
    TCGA_LUSC_PrimaryTumor_vs_TissueNormal <- read.csv("./output/TCGA_LUSC_PrimaryTumor_vs_TissueNormal_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    TCGA_LUSC_PrimaryTumor_vs_TissueNormal$dataset <- "TCGA_LUSC_PrimaryTumor_vs_TissueNormal"
    
    TCGA_READ_PrimaryTumor_vs_TissueNormal <- read.csv("./output/TCGA_READ_PrimaryTumor_vs_TissueNormal_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    TCGA_READ_PrimaryTumor_vs_TissueNormal$dataset <- "TCGA_READ_PrimaryTumor_vs_TissueNormal"
    
    TCGA_THCA_PrimaryTumor_vs_TissueNormal <- read.csv("./output/TCGA_THCA_PrimaryTumor_vs_TissueNormal_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    TCGA_THCA_PrimaryTumor_vs_TissueNormal$dataset <- "TCGA_THCA_PrimaryTumor_vs_TissueNormal"
    
    TCGA_BRCA_paired_PrimaryTumor_vs_TissueNormal <- read.csv("./output/TCGA_BRCA_paired_PrimaryTumor_vs_TissueNormal_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
    TCGA_BRCA_paired_PrimaryTumor_vs_TissueNormal$dataset <- "TCGA_BRCA_paired_PrimaryTumor_vs_TissueNormal"
    
    forbubbleplot <- rbind(GSE142978_HCC_vs_NC,GSE174302_CRC_vs_NC,GSE174302_ESCA_vs_NC,
                           GSE174302_HCC_vs_NC,GSE174302_LUAD_vs_NC,GSE174302_STAD_vs_NC,
                           GSE68086_BRCA_vs_NC,GSE68086_CRC_vs_NC,GSE68086_GBM_vs_NC,
                           GSE68086_Hepatobiliary_vs_NC,GSE68086_NSCLC_vs_NC,GSE68086_PDAC_vs_NC,
                           GSE89843_Multiple_Sclerosis_vs_NC,GSE89843_Pulmonary_Hypertension_vs_NC,GSE89843_Chronic_Pancreatitis_vs_NC,
                           GSE89843_Epilepsy_vs_NC,GSE89843_Atherosclerosis_vs_NC,GSE89843_Angina_Pectoris_vs_NC,GSE89843_NSCLC_vs_NC,
                           exoRBase_HCC_vs_NC,exoRBase_CRC_vs_NC,exoRBase_PAAD_vs_NC,
                           GSE133684_PDAC_vs_NC,GSE40174_PDAC_vs_NC,
                           GSE120663_HCC_vs_NC,GSE117623_HCC_vs_DC,
                           TCGA_COAD_PrimaryTumor_vs_TissueNormal,TCGA_STAD_PrimaryTumor_vs_TissueNormal,TCGA_BRCA_PrimaryTumor_vs_TissueNormal,
                           TCGA_ESCA_PrimaryTumor_vs_TissueNormal,TCGA_LIHC_PrimaryTumor_vs_TissueNormal,TCGA_LUAD_PrimaryTumor_vs_TissueNormal,
                           TCGA_LUSC_PrimaryTumor_vs_TissueNormal,TCGA_READ_PrimaryTumor_vs_TissueNormal,TCGA_THCA_PrimaryTumor_vs_TissueNormal)
    
    GSEA_res <- TCGA_COAD_PrimaryTumor_vs_TissueNormal
    filtered_up <- GSEA_res[(GSEA_res$p.adjust < 0.1) & (GSEA_res$NES > 0), ]
    up <- head(filtered_up[order(filtered_up$NES,decreasing = TRUE),]$Description,50)
    up <- factor(up,levels=up)
    filtered_down <- GSEA_res[(GSEA_res$p.adjust < 0.1) & (GSEA_res$NES < 0), ]
    down <- head(filtered_down[order(filtered_down$NES,decreasing = FALSE),]$Description,20)
    down <- factor(down,levels=down)
    top20 <- fct_c(up,down)
    
    genes <- top20
    genes <- c(
      "Metabolics","Glycolysis Gluconeogenesis","Lipoprotein metabolism","Biosynthesis of amino acids","Protein export","Protein secretion","Protein digestion absorption",
      "Pro-inflamatory cytokines&chemokines","Inflammatory response","Unfolded protein response","Primary immunodeficiency",
      "Epithelial mesenchymal transition","G2M checkpoint","Reactive oxigen species","Hypoxia","Angiogensis",
      "Apoptosis","Resisting Cell Death","Sustaining Proliferative Signaling","Extension of telomeres",
      "MT-RP-RNA","RP-mRNA","Ribosome"
    )
    
    
    forbubbleplot_plot <- filter(forbubbleplot, Description %in% genes)
    forbubbleplot_plot[forbubbleplot_plot$p.adjust>=0.05,"p.adjust"] <- NA
    
    forbubbleplot_plot$Description <- factor(forbubbleplot_plot$Description,levels = genes)
    forbubbleplot_plot$dataset <- factor(forbubbleplot_plot$dataset,levels = c("TCGA_COAD_PrimaryTumor_vs_TissueNormal","TCGA_STAD_PrimaryTumor_vs_TissueNormal","TCGA_BRCA_PrimaryTumor_vs_TissueNormal",
                                                                               "TCGA_ESCA_PrimaryTumor_vs_TissueNormal","TCGA_LIHC_PrimaryTumor_vs_TissueNormal","TCGA_LUAD_PrimaryTumor_vs_TissueNormal",
                                                                               "TCGA_LUSC_PrimaryTumor_vs_TissueNormal","TCGA_READ_PrimaryTumor_vs_TissueNormal","TCGA_THCA_PrimaryTumor_vs_TissueNormal",
                                                                               "GSE40174_PDAC_vs_NC","GSE120663_HCC_vs_NC","GSE117623_HCC_vs_DC",
                                                                               "exoRBase_HCC_vs_NC","exoRBase_CRC_vs_NC","exoRBase_PAAD_vs_NC","GSE133684_PDAC_vs_NC",
                                                                               "GSE142978_HCC_vs_NC","GSE174302_CRC_vs_NC","GSE174302_ESCA_vs_NC",
                                                                               "GSE174302_HCC_vs_NC","GSE174302_LUAD_vs_NC","GSE174302_STAD_vs_NC",
                                                                               "GSE68086_BRCA_vs_NC","GSE68086_CRC_vs_NC","GSE68086_GBM_vs_NC",
                                                                               "GSE68086_Hepatobiliary_vs_NC","GSE68086_NSCLC_vs_NC","GSE68086_PDAC_vs_NC","GSE89843_NSCLC_vs_NC",
                                                                               "GSE89843_Multiple-Sclerosis_vs_NC","GSE89843_Pulmonary-Hypertension_vs_NC","GSE89843_Chronic-Pancreatitis_vs_NC",
                                                                               "GSE89843_Epilepsy_vs_NC","GSE89843_Atherosclerosis_vs_NC","GSE89843_Angina-Pectoris_vs_NC"
                                                                               ))
    ggplot(forbubbleplot_plot,aes(x=forbubbleplot_plot$dataset,y=forbubbleplot_plot$Description))+
      geom_point(aes(size=-1*log10(forbubbleplot_plot$p.adjust),color=forbubbleplot_plot$NES))+
      scale_colour_gradient2(low="#162252",high="#8B0000",mid="white",midpoint = 0)+
      geom_vline(xintercept = c(9.5,12.5,16.5,22.5),linetype = "dashed",color="black")+
      geom_vline(xintercept = 29.5,linetype = "dashed",color="grey")+
      labs(color="NES",
           size=expression(-log[10](p.adjust)),
           x="Dataset")+
      theme_bw()+
      theme(axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
            axis.text.x = element_text(size=rel(1.3),face="bold",colour = "black", angle = 90, hjust = 1,vjust=0.5),
            axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
            axis.title.y = element_blank())
  }
}


#Figure1 sample cohort
{
##multi-omics sample shaozhen
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/01.Sample summary")   

#sunburst
p <- read.csv("Figure3_forplot.csv",header = F, stringsAsFactors = FALSE)
c <- read.csv("Color.csv",header = F,stringsAsFactors=FALSE)
c[1,1] <- "NA"
sunburst(p,
         count=TRUE,
         sortFunction = htmlwidgets::JS(
           "function(a,b) {
                       // sort by count descending
                       // unlike the other example using data.name, value is at the top level of the object
                       return b.value - a.value}" ),
         legend=list(w=240,h=30,r=10,s=5),
         #sortFunction = htmlwidgets::JS("Health","Colorectum","Liver","Stomach","Lung","Thyroid","Esophagus"),
         legendOrder = list("NC","CRC","STAD","Female","Male","Stage I","Stage II","Stage III","Stage IV","Unknown","long RNA","DNA methylation","WGS","small RNA"),
         colors=list(range = c$V2, domain = c$V1)
)  

#Age_Stage barplot
clinical <- read.csv("./Figure3/Figure3_subtypes_forplot.csv",header = TRUE)
View(clinical)
clinical[clinical==""]<-NA

clinical$Stage.summary <- factor(clinical$Stage.summary,levels = c("Stage I","Stage II","Stage III","Stage IV"))
ggplot(data=clinical[which(clinical$Type!="NC"),],aes(x=Age_10,y=number,fill=Stage.summary))+ geom_bar(stat = "identity", width=0.9, col='transparent') +
  theme_bw()+
  theme(#legend.position="bottom",
    legend.position="right",
    panel.grid.major.x = element_blank(),
    legend.title = element_text(face="bold", color="transparent",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(color="black", size=20,angle = 0,hjust = 0.5),
    axis.text.y = element_text(color="black", size=20),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  ylab("Number of patients")+xlab("Ages(years)")+
  scale_y_continuous(expand = c(0,0),limits = c(0,35))+
  scale_fill_manual(values=c("#DCDCDC","#A3A3A3","#4A4A4A","#000000"))



#subtype
clinical <- read.csv("./Figure3/Figure3_simple_subtype_forplot.csv",header = TRUE,stringsAsFactors = FALSE)
color <- read.csv("Color_for_subtype.csv",header = FALSE)
clinical[clinical==""]<-NA
clinical[is.na(clinical)] <- "No biopsy"

color$V1 <- factor(color$V1,levels = color$V1)

clinical$index <- factor(clinical$index,levels = c("CRC position","MMR","STAD position","HER2"))
clinical$Subtype <- factor(clinical$Subtype,levels = color$V1)
ggplot(data=clinical,aes(x=number,y=number,fill=clinical$Subtype))+ geom_bar(stat = "identity", width=10, col='transparent') + facet_grid(~index)+
  #geom_text(aes(label=Subtype), position = "stack",vjust=1)+
  theme_bw()+
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.title = element_text(face="bold", color="transparent",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(color="black", size=20,angle = 0,hjust = 0.5),
    axis.text.y = element_text(color="black", size=20),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24),
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=10))+
  ylab("Number of patients")+xlab("")+scale_y_reverse()+scale_x_discrete(position = "top")+
  scale_fill_manual(values=c("#FFFFFF","#DCDCDC","#FFF8DC","#EEDD82","#EDCB62","#FCB514","#FEF1B5","#E5BC3B","#CD9B1D","#EDCB62","#FCB514","#E5BC3B","#CD9B1D"))
  #scale_fill_manual(values=c("#FFFFFF","#DCDCDC","#FFF8DC","#EEDD82","#FEE5AC","#EDCB62","#FCB514","#FEF1B5","#E5BC3B","#CD9B1D","#EDCB62","#FCB514","#E5BC3B","#CD9B1D"))
}

#Sample QC
{
  library(ggplot2)
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/quality control")
  qc <- read.csv("quality_control_forplot.csv",header = TRUE)
  
  qc$Criterian <- factor(qc$Criterian,level=c("Criterian 8: Intron spanning","Criterian 7: Unclassified ratio","Criterian 6: Long RNA ratio","Criterian 5: rRNA ratio","Criterian 4: Genome aligned","Criterian 3: Spike in ratio","Criterian 2: Clean","Criterian 1: Raw"))
  qc$Results <- factor(qc$Results,level=c("TRUE","FALSE"))
  ggplot(qc,aes(x=Criterian,y=log10(qc$X.Read.pairs.or.Ratio..Criterien)))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(0.5),outlier.size=0,outlier.alpha = 0)+
    geom_point(aes(color = qc$Results),size = 1, position = position_jitterdodge(dodge.width=0.5,jitter.width = 0.5))+
    scale_color_manual(values=c("black","red")) +
    theme_bw()+
    coord_flip()+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(size=1, colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20,angle = 0,vjust = 1, hjust=0.5),
      axis.text.y = element_text(face="bold",  color="black", size=18,vjust = 0.5, hjust=0),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    geom_hline(aes(yintercept=0),linetype=5,col="grey")+
    labs(color="Passed")+
    ylab("Log10((Read pairs or Ratio)/Criterien)")+
    xlab("")
}

#normalized gene coverage between RNA and DNA
{
{
  {
  library(ggbreak) 
  library(patchwork)
  target_gene <- "APC"
  setwd(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Coverage/CRC-PKU-32_normalized/",target_gene))
  DNA <- read.csv("CRC-PKU-32-wgs_sorted_cpm.txt",sep = "\t",header = FALSE)
  methy <- read.csv("CRC-PKU-32-me_sorted_cpm.txt",sep = "\t",header = FALSE)
  #RNA <- read.csv("CRC-PKU-32-pico-intron-spanning_sorted_cpm.txt",sep = "\t",header = FALSE)
  RNA <- read.csv("CRC-PKU-32-genome_rmdup_sorted_cpm.txt",sep = "\t",header = FALSE)
  
  if(file.info("CRC-PKU-32-qia_sorted_cpm.txt")$size == 0) {
    qia <- as.data.frame(matrix(numeric(0),ncol=3))
  } else {
    qia <- read.csv("CRC-PKU-32-qia_sorted_cpm.txt",sep = "\t",header = FALSE)
  }
  
  colnames(DNA) <- c("Chr","Position","Depth")
  colnames(methy) <- c("Chr","Position","Depth")
  colnames(RNA) <- c("Chr","Position","Depth")
  colnames(qia) <- c("Chr","Position","Depth")
  
  RNA_DNA <- full_join(RNA,DNA,by = c("Position"="Position"))
  
  colnames(RNA_DNA) <- c("Chr.RNA","Position","Depth.RNA","Chr.DNA","Depth.DNA") 
  
  RNA_DNA_methy <- full_join(RNA_DNA,methy,by = c("Position"="Position"))
  
  colnames(RNA_DNA_methy) <- c("Chr.RNA","Position","Depth.RNA","Chr.DNA","Depth.DNA","Chr.methy","Depth.methy")
  
  RNA_DNA_methy_qia <- full_join(RNA_DNA_methy,qia,by = c("Position"="Position"))
  
  colnames(RNA_DNA_methy_qia) <- c("Chr.RNA","Position","Depth.RNA","Chr.DNA","Depth.DNA","Chr.methy","Depth.methy","Chr.qia","Depth.aiq")
  
  RNA_DNA_methy_qia <- RNA_DNA_methy_qia[,-8]
  RNA_DNA_methy_qia <- RNA_DNA_methy_qia[,-6]
  RNA_DNA_methy_qia <- RNA_DNA_methy_qia[,-4]
  RNA_DNA_methy_qia <- RNA_DNA_methy_qia[,-1]
  
  RNA_DNA_methy_qia[is.na(RNA_DNA_methy_qia)] <- 0
  
  RNA_forplot <- RNA_DNA_methy_qia[,-c(3,4,5)]
  RNA_forplot$type <- "Total cfRNA"
  colnames(RNA_forplot) <- c("Position","Depth","type")
  
  DNA_forplot <- RNA_DNA_methy_qia[,-c(2,4,5)]
  DNA_forplot$type <- "cfDNA WGS"
  colnames(DNA_forplot) <- c("Position","Depth","type")
  
  
  methy_forplot <- RNA_DNA_methy_qia[,-c(2,3,5)]
  methy_forplot$type <- "cfDNA methylation"
  colnames(methy_forplot) <- c("Position","Depth","type")
  
  qia_forplot <- RNA_DNA_methy_qia[,-c(2,3,4)]
  qia_forplot$type <- "Small cfRNA"
  colnames(qia_forplot) <- c("Position","Depth","type")
  
  
  
  plot <- rbind(RNA_forplot,DNA_forplot,methy_forplot,qia_forplot)
  
  gene_loc <- read.csv(paste0(target_gene,".bed"),sep="\t",header=FALSE)
  up_stream <- gene_loc[,2]
  down_stream <- gene_loc[,3]
  
  plot <- plot[plot$Position>up_stream,]
  plot <- plot[plot$Position<down_stream,]
  plot$layer2 <- plot$Depth
  plot[-grep("Total cfRNA",plot$type),]$layer2 <- NA
  plot$layer3 <- plot$Depth
  plot[-grep("Small cfRNA",plot$type),]$layer3 <- NA
  plot$layer4 <- plot$Depth
  plot[-grep("cfDNA methylation",plot$type),]$layer4 <- NA

  plot$type <- factor(plot$type,levels=c("Total cfRNA","cfDNA WGS","cfDNA methylation","Small cfRNA"))
  plot <- plot[(plot$Position>112834750)&(plot$Position<112846750),]
  
  RNA_upperlimit <- ceiling(max(plot[grep("Total cfRNA",plot$type),]$Depth)/0.1)*0.1
  WGS_meanCoverage <- ceiling(mean(plot[grep("cfDNA WGS",plot$type),]$Depth)/0.1)*0.1
  #WGS_meanCoverage <- 20
  
  p <- ggplot(plot,aes(x=plot$Position,y=log(plot$Depth+1,10),fill=plot$type))+
    geom_bar(stat="identity",position = 'dodge')+
    #geom_bar(aes(x=plot$Position,y=log10(plot$layer2+1),fill=plot$type),stat="identity",position = 'dodge')+
    geom_bar(aes(x=plot$Position,y=log10(plot$layer3+1),fill=plot$type),stat="identity",position = 'dodge')+
    geom_bar(aes(x=plot$Position,y=log10(plot$layer4+1),fill=plot$type),stat="identity",position = 'dodge')+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(limits=c(0, log(1+RNA_upperlimit+RNA_upperlimit*0.8,10)),
                       breaks = c(0,log(WGS_meanCoverage+1,10),log(RNA_upperlimit+1,10)),
                       labels=c("",WGS_meanCoverage,RNA_upperlimit),expand = c(0,0))+
    #scale_x_discrete(expand=c(0.2, 0.2))+
    #scale_fill_jco(alpha = 0.8)+
    #scale_fill_manual(values=c(alpha("#5CACEE","#FFCC11",alpha = 0.45),"#666666","#CD0000"))+
    scale_fill_manual(values=c("#5CACEE",alpha("#EEB4B4",alpha = 0.45),"#666666","#CD0000"))+
    theme_bw()+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(size=1, colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=18,angle = 0,hjust = 0.5,vjust=0.5),
      axis.text.y = element_text(face="bold",  color="black", size=24),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=23))+
    xlab("")+
    ylab("(Base depth)/(million read pairs))")+
    geom_hline(yintercept=log(WGS_meanCoverage+1,10),color = alpha("red",alpha = 0.7),size = 1.5,linetype = 2)+
    labs(fill="Type")
  
  #p <- p+scale_y_break(c(log(0.5+1,10), log(1+RNA_upperlimit-RNA_upperlimit*0.8,10)),
  #                     scale = 0.7,ticklabels = NULL)+
  #  scale_y_continuous(limits=c(0, log(1+RNA_upperlimit+RNA_upperlimit*0.8,10)),
  #                     breaks = c(0,log(WGS_meanCoverage+1,10),log(RNA_upperlimit-100+1,10)),
  #                     labels=c("",WGS_meanCoverage,RNA_upperlimit-100),expand = c(0,0))
  
  ggsave(paste0(target_gene,"_genome_rmdup.pdf"),p,width=18,height=6)
  }
  
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(trackViewer)
  #GAPDH: chr12	6533926	6538374	ENSG00000111640.14	0	+
  #APC: chr5	112707497	112846239	ENSG00000134982.16	0	+
  #KRAS: chr12	25204788	25250936	ENSG00000133703.11	0	-
  #TP53: chr17	7661778	7687550	ENSG00000141510.16	0	-
  gr <- GRanges("chr5", IRanges(112834750, 112846750), strand="+")
  trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                           org.Hs.eg.db,
                           gr=gr)
  library(trackViewer)
  viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
  trackList <- trackList(trs)
  vp <- viewTracks(trackList, 
                   gr=gr, viewerStyle=viewerStyle, 
                   autoOptimizeStyle=TRUE)
  addGuideLine(c(6533927, 6538374), vp=vp)
  #addArrowMark(list(x=6533927, 
  #                  y=1), # 2 means track 2 from the bottom.
  #             label="label",
  #             col="dark blue",
  #             vp=vp)
}

{
  library(Gviz)
  data(geneModels)
  head(geneModels)
  nrow(geneModels)
  
  Track <- BiomartGeneRegionTrack(chromosome = 7, 
                                  start=6534434, end=6538374, name = "ENSEMBL", 
                                  biomart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl"))
  Track <- BiomartGeneRegionTrack(
    start=6534434, end=6538374,chromosome="chr7",genome='hg38')
  
  library(rtracklayer)
  library(trackViewer)
  extdata <- system.file("extdata", package="trackViewer",
                         mustWork=TRUE)
  gr <- GRanges("chr12", IRanges(6533927, 6538374), strand="-")
  GAPDH <- importScore(file.path("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/Coverage/GAPDH.bed"), format="BED",
                       ranges=gr)
  #fox2$dat <- coverageGR(fox2$dat)
  
  viewTracks(trackList(GAPDH), gr=gr, autoOptimizeStyle=TRUE, newpage=FALSE)
}
}

#feature summary
{
  library(Rmisc)
  library(dplyr)
  library(ggplot2)
  library(ggsci)
  feature_summary={}
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features")
  Alt.promoter <- read.csv("./Alt.promoter/TPM-by-tx_multiomics_paired.txt",sep="\t",header = TRUE,row.names=1)
  Alt.promoter[is.na(Alt.promoter)] <- 0
  feature_Alt.promoter <- as.data.frame(colSums(Alt.promoter>0))
  colnames(feature_Alt.promoter) <- "Number"
  feature_Alt.promoter$feature <- "Alt.promoter"
  feature_Alt.promoter$sample <- rownames(feature_Alt.promoter)
  feature_Alt.promoter$group <- as.character(lapply(strsplit(feature_Alt.promoter$sample,".",fixed = TRUE),function(x) x[1]))
  head(feature_Alt.promoter)
  feature_summary <- rbind(feature_summary,feature_Alt.promoter)
  #median(feature_Alt.promoter[grep("NC",rownames(feature_Alt.promoter)),])
  #median(feature_Alt.promoter[grep("CRC",rownames(feature_Alt.promoter)),])
  #median(feature_Alt.promoter[grep("STAD",rownames(feature_Alt.promoter)),])
  
  APA <- read.csv("./APA/PDUI.txt",sep="\t",header = TRUE,row.names=1)
  APA[is.na(APA)] <- 0
  feature_APA <- as.data.frame(colSums(APA>0))
  colnames(feature_APA) <- "Number"
  feature_APA$feature <- "APA"
  feature_APA$sample <- rownames(feature_APA)
  feature_APA$group <- as.character(lapply(strsplit(feature_APA$sample,".",fixed = TRUE),function(x) x[1]))
  head(feature_APA)
  feature_summary <- rbind(feature_summary,feature_APA)
  #median(feature_APA[grep("NC",rownames(feature_APA)),])
  #sd(feature_APA[grep("NC",rownames(feature_APA)),])
  #median(feature_APA[grep("CRC",rownames(feature_APA)),])
  #sd(feature_APA[grep("CRC",rownames(feature_APA)),])
  #median(feature_APA[grep("STAD",rownames(feature_APA)),])
  #sd(feature_APA[grep("STAD",rownames(feature_APA)),])
  
  Chimeric <- read.csv("./Chimeric RNA/multiomics_paired_chimeric.txt",sep="\t",header = TRUE,row.names=1)
  Chimeric[is.na(Chimeric)] <- 0
  feature_Chimeric <- as.data.frame(colSums(Chimeric>0))
  colnames(feature_Chimeric) <- "Number"
  feature_Chimeric$feature <- "Chimeric"
  feature_Chimeric$sample <- rownames(feature_Chimeric)
  feature_Chimeric$group <- as.character(lapply(strsplit(feature_Chimeric$sample,".",fixed = TRUE),function(x) x[1]))
  head(feature_Chimeric)
  feature_summary <- rbind(feature_summary,feature_Chimeric)
  #median(feature_Chimeric[grep("NC",rownames(feature_Chimeric)),])
  #median(feature_Chimeric[grep("CRC",rownames(feature_Chimeric)),])
  #median(feature_Chimeric[grep("STAD",rownames(feature_Chimeric)),])
  
  Microbe <- read.csv("./microbe/microbe_G_level.txt",sep="\t",header = TRUE,row.names=1)
  Microbe[is.na(Microbe)] <- 0
  feature_Microbe <- as.data.frame(colSums(Microbe>0))
  colnames(feature_Microbe) <- "Number"
  feature_Microbe$feature <- "Microbe_Genus"
  feature_Microbe$sample <- rownames(feature_Microbe)
  feature_Microbe$group <- as.character(lapply(strsplit(feature_Microbe$sample,".",fixed = TRUE),function(x) x[1]))
  head(feature_Microbe)
  feature_summary <- rbind(feature_summary,feature_Microbe)
  #median(feature_Microbe[grep("NC",rownames(feature_Microbe)),])
  #median(feature_Microbe[grep("CRC",rownames(feature_Microbe)),])
  #median(feature_Microbe[grep("STAD",rownames(feature_Microbe)),])
  
  TE <- read.csv("./TE/multiomics_paired_TE_readnum.txt",sep="\t",header = TRUE,row.names=1)
  TE <- TE[-1,]
  TE2 <- mutate_all(TE, function(x) as.numeric(as.character(x)))
  rownames(TE2) <- rownames(TE)
  TE2[is.na(TE2)] <- 0
  feature_TE <- as.data.frame(colSums(TE2>0))
  colnames(feature_TE) <- "Number"
  feature_TE$feature <- "TE"
  feature_TE$sample <- rownames(feature_TE)
  feature_TE$group <- as.character(lapply(strsplit(feature_TE$sample,".",fixed = TRUE),function(x) x[1]))
  head(feature_TE)
  feature_summary <- rbind(feature_summary,feature_TE)
  #median(feature_TE[grep("NC",rownames(feature_TE)),])
  #median(feature_TE[grep("CRC",rownames(feature_TE)),])
  #median(feature_TE[grep("STAD",rownames(feature_TE)),])
  
  Splicing_A3SS <- read.csv("./splicing/A3SS_JC_inc.txt",sep="\t",header = TRUE,row.names=1)
  colnames(Splicing_A3SS) <- gsub("X.data.taoyuhuan.projects.exOmics_RNA.multiomics_paired.output.Intron.spanning.","",colnames(Splicing_A3SS))
  colnames(Splicing_A3SS) <- gsub(".intron.spanning.bam","",colnames(Splicing_A3SS))
  Splicing_A3SS[is.na(Splicing_A3SS)] <- 0
  feature_Splicing_A3SS <- as.data.frame(colSums(Splicing_A3SS>0))
  colnames(feature_Splicing_A3SS) <- "Number"
  feature_Splicing_A3SS$feature <- "Splicing_A3SS"
  feature_Splicing_A3SS$sample <- rownames(feature_Splicing_A3SS)
  feature_Splicing_A3SS$group <- as.character(lapply(strsplit(feature_Splicing_A3SS$sample,".",fixed = TRUE),function(x) x[1]))
  head(feature_Splicing_A3SS)
  feature_summary <- rbind(feature_summary,feature_Splicing_A3SS)
  #median(feature_Splicing_A3SS[grep("NC",rownames(feature_Splicing_A3SS)),])
  #median(feature_Splicing_A3SS[grep("CRC",rownames(feature_Splicing_A3SS)),])
  #median(feature_Splicing_A3SS[grep("STAD",rownames(feature_Splicing_A3SS)),])
  
  Splicing_A5SS <- read.csv("./splicing/A5SS_JC_inc.txt",sep="\t",header = TRUE,row.names=1)
  colnames(Splicing_A5SS) <- gsub("X.data.taoyuhuan.projects.exOmics_RNA.multiomics_paired.output.Intron.spanning.","",colnames(Splicing_A5SS))
  colnames(Splicing_A5SS) <- gsub(".intron.spanning.bam","",colnames(Splicing_A5SS))
  Splicing_A5SS[is.na(Splicing_A5SS)] <- 0
  feature_Splicing_A5SS <- as.data.frame(colSums(Splicing_A5SS>0))
  colnames(feature_Splicing_A5SS) <- "Number"
  feature_Splicing_A5SS$feature <- "Splicing_A5SS"
  feature_Splicing_A5SS$sample <- rownames(feature_Splicing_A5SS)
  feature_Splicing_A5SS$group <- as.character(lapply(strsplit(feature_Splicing_A5SS$sample,".",fixed = TRUE),function(x) x[1]))
  head(feature_Splicing_A5SS)
  feature_summary <- rbind(feature_summary,feature_Splicing_A5SS)
  #median(feature_Splicing_A5SS[grep("NC",rownames(feature_Splicing_A5SS)),])
  #median(feature_Splicing_A5SS[grep("CRC",rownames(feature_Splicing_A5SS)),])
  #median(feature_Splicing_A5SS[grep("STAD",rownames(feature_Splicing_A5SS)),])
  
  Splicing_MXE <- read.csv("./splicing/MXE_JC_inc.txt",sep="\t",header = TRUE,row.names=1)
  colnames(Splicing_MXE) <- gsub("X.data.taoyuhuan.projects.exOmics_RNA.multiomics_paired.output.Intron.spanning.","",colnames(Splicing_MXE))
  colnames(Splicing_MXE) <- gsub(".intron.spanning.bam","",colnames(Splicing_MXE))
  Splicing_MXE[is.na(Splicing_MXE)] <- 0
  feature_Splicing_MXE <- as.data.frame(colSums(Splicing_MXE>0))
  colnames(feature_Splicing_MXE) <- "Number"
  feature_Splicing_MXE$feature <- "Splicing_MXE"
  feature_Splicing_MXE$sample <- rownames(feature_Splicing_MXE)
  feature_Splicing_MXE$group <- as.character(lapply(strsplit(feature_Splicing_MXE$sample,".",fixed = TRUE),function(x) x[1]))
  head(feature_Splicing_MXE)
  feature_summary <- rbind(feature_summary,feature_Splicing_MXE)
  #median(feature_Splicing_MXE[grep("NC",rownames(feature_Splicing_MXE)),])
  #median(feature_Splicing_MXE[grep("CRC",rownames(feature_Splicing_MXE)),])
  #median(feature_Splicing_MXE[grep("STAD",rownames(feature_Splicing_MXE)),])
  
  Splicing_RI <- read.csv("./splicing/RI_JC_inc.txt",sep="\t",header = TRUE,row.names=1)
  colnames(Splicing_RI) <- gsub("X.data.taoyuhuan.projects.exOmics_RNA.multiomics_paired.output.Intron.spanning.","",colnames(Splicing_RI))
  colnames(Splicing_RI) <- gsub(".intron.spanning.bam","",colnames(Splicing_RI))
  Splicing_RI[is.na(Splicing_RI)] <- 0
  feature_Splicing_RI <- as.data.frame(colSums(Splicing_RI>0))
  colnames(feature_Splicing_RI) <- "Number"
  feature_Splicing_RI$feature <- "Splicing_RI"
  feature_Splicing_RI$sample <- rownames(feature_Splicing_RI)
  feature_Splicing_RI$group <- as.character(lapply(strsplit(feature_Splicing_RI$sample,".",fixed = TRUE),function(x) x[1]))
  head(feature_Splicing_RI)
  feature_summary <- rbind(feature_summary,feature_Splicing_RI)
  #median(feature_Splicing_RI[grep("NC",rownames(feature_Splicing_RI)),])
  #median(feature_Splicing_RI[grep("CRC",rownames(feature_Splicing_RI)),])
  #median(feature_Splicing_RI[grep("STAD",rownames(feature_Splicing_RI)),])
  
  Splicing_SE <- read.csv("./splicing/SE_JC_inc.txt",sep="\t",header = TRUE,row.names=1)
  colnames(Splicing_SE) <- gsub("X.data.taoyuhuan.projects.exOmics_RNA.multiomics_paired.output.Intron.spanning.","",colnames(Splicing_SE))
  colnames(Splicing_SE) <- gsub(".intron.spanning.bam","",colnames(Splicing_SE))
  Splicing_SE[is.na(Splicing_SE)] <- 0
  feature_Splicing_SE <- as.data.frame(colSums(Splicing_SE>0))
  colnames(feature_Splicing_SE) <- "Number"
  feature_Splicing_SE$feature <- "Splicing_SE"
  feature_Splicing_SE$sample <- rownames(feature_Splicing_SE)
  feature_Splicing_SE$group <- as.character(lapply(strsplit(feature_Splicing_SE$sample,".",fixed = TRUE),function(x) x[1]))
  head(feature_Splicing_SE)
  feature_summary <- rbind(feature_summary,feature_Splicing_SE)
  #median(feature_Splicing_SE[grep("NC",rownames(feature_Splicing_SE)),])
  #median(feature_Splicing_SE[grep("CRC",rownames(feature_Splicing_SE)),])
  #median(feature_Splicing_SE[grep("STAD",rownames(feature_Splicing_SE)),])
  
  
  #Splicing_A3SS_stat <- read.csv("./splicing/A3SS_JC_stats.txt",sep = "\t", header = TRUE, row.names = 1)
  #nrow(Splicing_A3SS_stat[Splicing_A3SS_stat$FDR<0.05,])
  
  #Splicing_A5SS_stat <- read.csv("./splicing/A5SS_JC_stats.txt",sep = "\t", header = TRUE, row.names = 1)
  #nrow(Splicing_A5SS_stat[Splicing_A5SS_stat$FDR<0.05,])
  
  #Splicing_MXE_stat <- read.csv("./splicing/MXE_JC_stats.txt",sep = "\t", header = TRUE, row.names = 1)
  #nrow(Splicing_MXE_stat[Splicing_MXE_stat$FDR<0.05,])
  
  #Splicing_RI_stat <- read.csv("./splicing/RI_JC_stats.txt",sep = "\t", header = TRUE, row.names = 1)
  #nrow(Splicing_RI_stat[Splicing_RI_stat$FDR<0.05,])
  
  #Splicing_SE_stat <- read.csv("./splicing/SE_JC_stats.txt",sep = "\t", header = TRUE, row.names = 1)
  #nrow(Splicing_SE_stat[Splicing_SE_stat$FDR<0.05,])
  
  #ASE and Editing are added manully
  write.csv(feature_summary,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/feature_summary.csv")
  feature_summary <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Old/feature_summary.csv",header=TRUE)
  
  feature_summary_count <- summarySE(feature_summary,measurevar = "Number",groupvars = c("feature","group"))
  

  feature_summary_count$group <- factor(feature_summary_count$group,levels = c("CRC","STAD","NC"))
  feature_summary_count$feature <- factor(feature_summary_count$feature,levels = levels(as.factor(feature_summary_count$feature)))
  ggplot(feature_summary_count, aes(x = feature, y = log10(Number),fill = feature)) +
    facet_wrap(~group)+
    coord_flip()+
    geom_bar(stat = "identity", color = "black", position = position_dodge()) +
    geom_errorbar(aes(ymin = log10(Number - ci), ymax = log10(Number + ci)), 
                  width = 0.2, position = position_dodge(0.9))+
    scale_fill_simpsons()+
    scale_y_continuous(breaks=c(0,1,2,3,4,5),labels=c("","10","","1000","","100000"))+
    theme_bw()+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(size=0.5, colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
      strip.text = element_text(face="bold",  color="black", size=24),
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_text(color="black", size=16,angle = 0,hjust = 0.5),
      axis.text.y = element_text(face="bold",  color="black", size=24),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
}

#diff feature summary
{
  Features_color <- c("DNA: copy number"="#FFE600",
                      "DNA: nucleosome occupancy"="#60AFFE",
                      "DNA: methylation region"="#CDC9C9",
                      "RNA: alternative promoter"="#E9C2A6",
                      "RNA: expression"="#A5435C",
                      "RNA: splicing event"="#C5E3BF",
                      "RNA: alternative polyadenyltion"="#003F87",
                      "RNA: chimeric event"="#FF3D0D",
                      "RNA: editing event"="#324F17",
                      "RNA: allele specific expression"="#87CEFF",
                      "RNA: SNP"="#333333")
  
  
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20210722_multiomics/")
  
  all_features <- c("DNA: copy number"="DNA-CNV",
                    "DNA: nucleosome occupancy"="DNA-nucleosome",
                    "DNA: methylation region"="MeDIP",
                    "RNA: alternative promoter"="Alt.promoter",
                    "RNA: expression"="Expression",
                    "RNA: SNP"="SNP",
                    "RNA: allele specific expression"="ASE",
                    "RNA: splicing event"="Splicing_FDR0.05",
                    "RNA: chimeric event"="chimeric",
                    "RNA: editing event"="EDIT",
                    "RNA: alternative polyadenyltion"="APA")
  
  i=1
  j=1
  Group <- "STADvsNC"
  diff_feature <- as.data.frame(matrix(numeric(0),ncol=4))
  colnames(diff_feature) <- c("Feature","Group","Number","Type")
  
  while (i <= length(all_features)){
  diff <- read.csv(paste0(all_features[i],"/",Group,"/Diff_wilcox.txt"),header = TRUE,row.names = 1,sep="\t")
  up <- nrow(diff[(diff$pvalue<0.001)&(diff$log2FoldChange>0),])
  down <- nrow(diff[(diff$pvalue<0.001)&(diff$log2FoldChange<0),])
  diff_feature[c(j,j+1),"Feature"] <- names(all_features)[i]
  diff_feature[c(j,j+1),"Group"] <- Group
  diff_feature[j,"Number"] <- up
  diff_feature[j,"Type"] <- "Up-regulated"
  diff_feature[j+1,"Number"] <- down
  diff_feature[j+1,"Type"] <- "Down-regulated"
  j=j+2
  i=i+1
  }
  
  CRC <- diff_feature
  STAD <- diff_feature
  diff_feature <- rbind(CRC,STAD)
  diff_feature <- aggregate(Number ~ Feature + Group, data = diff_feature, FUN = sum, na.rm = TRUE)
  
  diff_feature$Feature <- factor(diff_feature$Feature,levels = rev(names(all_features)))
  #diff_feature$Type <- factor(diff_feature$Type,levels = c("Down-regulated","Up-regulated"))
  ggplot(diff_feature, aes(x = Feature, y = log10(Number),fill = Feature)) +
    facet_wrap(~Group)+
    coord_flip()+
    geom_bar(stat = "identity", position = "stack", color = "black")+
    scale_fill_manual(values = Features_color)+
    #scale_color_manual(values = c(NA,"black"))+
    scale_y_continuous(breaks=c(log10(1),log10(10),log10(100),log10(4000)),labels=c("","10","100","4000"))+
    geom_vline(xintercept=8.5, linetype="dashed", color = "grey")+
    theme_bw()+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(size=0.5, colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
      strip.text = element_text(face="bold",  color="black", size=24),
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_text(color="black", size=16,angle = 0,hjust = 0.5),
      axis.text.y = element_blank(),
      #axis.text.y = element_text(face="bold",  color="black", size=16),
      axis.title.x = element_text(face="bold", color="black", size=16),
      axis.title.y = element_text(face="bold",color="black", size=16))
}

##waterfall plot by maftools
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("maftools")
  
  library(maftools)
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/07.SNPs/aggregation_specimens/")
  var.annovar.maf = annovarToMaf(annovar = "specimens_comparison_forplot.txt", 
                                 Center = 'NA', 
                                 refBuild = 'hg38', 
                                 tsbCol = 'Tumor_Sample_Barcode', 
                                 table = 'refGene',
                                 sep = "\t")
  write.table(var.annovar.maf,file="var_annovar.maf",quote= F,sep="\t",row.names=F)
  var_maf = read.maf(maf ="var_annovar.maf")
  
  var_maf@clinical.data$Type <- as.character(lapply(strsplit(as.character(var_maf@clinical.data$Tumor_Sample_Barcode),"-"),function(x) x[1]))
  
  var_maf@clinical.data
  
  plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median')
  
  type_colors = c("red","blue")
  names(type_colors) = c("CRC", "STAD")
  type_colors = list(FAB_classification = type_colors)
  
  oncoplot(maf = var_maf, top = 18,showTumorSampleBarcodes = F,clinicalFeatures = "Type",
           fontSize = 0.8,
           genes = c("Plasma-APC","Tumor-APC","PBMC-APC",
                     "Plasma-KRAS","Tumor-KRAS","PBMC-KRAS",
                     "Plasma-TP53","Tumor-TP53","PBMC-TP53",
                     "Plasma-B2M","Tumor-B2M","PBMC-B2M",
                     "Plasma-CUX1","Tumor-CUX1","PBMC-CUX1",
                     "Plasma-UBR5","Tumor-UBR5","PBMC-UBR5"),
           sampleOrder = c("CRC-PKU-27","CRC-PKU-28","CRC-PKU-29","CRC-PKU-30","CRC-PKU-31","CRC-PKU-32","CRC-PKU-34",
                           "CRC-PKU-35","CRC-PKU-36","CRC-PKU-37","CRC-PKU-38","CRC-PKU-39","CRC-PKU-40","CRC-PKU-41"),
           #genes = c("RNA-APC","DNA-APC","RNA-KRAS","DNA-KRAS","RNA-TP53","DNA-TP53","RNA-B2M","DNA-B2M","RNA-CUX1","DNA-CUX1","RNA-UBR5","DNA-UBR5"),
           #annotationColor = type_colors,
           #annotationOrder = c("cfRNA","PBMC","Tumor"),
           #annotationOrder = c("PAAD"),
           colbar_pathway = FALSE,
           sortByAnnotation = TRUE,
           keepGeneOrder = TRUE,
           draw_titv = TRUE)
  
  laml.titv = titv(maf = var_maf, plot = FALSE, useSyn = TRUE)
  plotTiTv(res = laml.titv)
  
  somaticInteractions(maf = var_maf, top = 10, pvalue = c(0.05, 0.1))
  
  #lollipop plot for APC
  lollipopPlot(
    maf = var_maf,
    gene = 'APC',
    AACol = 'aaChange',
    showMutationRate = TRUE,
    #labelPos = "all",
    #refSeqID = "NM_000546",
    printCount = TRUE,
    showDomainLabel = TRUE
  )
  
  #oncogene pathway
  OncogenicPathways(maf = var_maf)
  PlotOncogenicPathways(maf = var_maf,pathways = "WNT")
  
  #
  rainfallPlot(maf = var_maf, detectChangePoints = TRUE, pointSize = 0.4)
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
  
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/PCAWG_20210722/")
  candidate_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Colorectal Cancer Genes.csv",header = TRUE)
  #choose a confidence level
  level <- 0.05 #c(0.05,0.01)
  #choose an interested cancer type
  group <- "STAD" #c("CRC","STAD","panCancer")
  #choose a cutoff, then subset
  cutoff <- "BP_STAD_0.05" #c("BP_CRC_0.05","BP_CRC_0.01","BP_STAD_0.05","BP_STAD_0.01","BP_panCancer_0.05","BP_panCancer_0.01","CRC_up0.05","CRC_up0.01","CRC_down0.05","CRC_down0.01","STAD_up0.05","STAD_up0.01","STAD_down0.05","STAD_down0.01","panCancer_up0.05","panCancer_up0.01","panCancer_down0.05","panCancer_down0.01")
  
  #give the total sample number
  
  #CRC sample number
  #{
  #sampleN_RNA <- 23
  #sampleN_DNA <- 23
  #}
  #STAD sample number
  {
    sampleN_RNA <- 30
    sampleN_DNA <- 30
  }
  
  #panCancer sample number
  #{
  #  sampleN_RNA <- 61
  #  sampleN_DNA <- 61
  #}
  
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
  APA <- read.csv("APA/all.csv",header = TRUE)
  APA <- APA[,-1]
  rownames(APA) <- APA$gene
  APA <- select_bestperformance(APA,level,group)
  APA_targetgenes <- subset_row(APA,candidate_list,cutoff)
  
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
  #summarize outliers (RNA)
  {
    #input features
    Alt_cutoff <- Alt_targetgenes[grep(cutoff,colnames(Alt_targetgenes))]
    Expression_cutoff <- Expression_targetgenes[grep(cutoff,colnames(Expression_targetgenes))]
    Splicing_cutoff <- Splicing_targetgenes[grep(cutoff,colnames(Splicing_targetgenes))]
    
    APA_cutoff <- APA_targetgenes[grep(cutoff,colnames(APA_targetgenes))]
    chimeric_cutoff <- chimeric_targetgenes[grep(cutoff,colnames(chimeric_targetgenes))]
    
    Editing_cutoff <- Editing_targetgenes[grep(cutoff,colnames(Editing_targetgenes))]
    ASE_cutoff <- ASE_targetgenes[grep(cutoff,colnames(ASE_targetgenes))]
    SNP_cutoff <- SNP_targetgenes[grep(cutoff,colnames(SNP_targetgenes))]
    
    features <- cbind(Alt_cutoff,Expression_cutoff,Splicing_cutoff,APA_cutoff,chimeric_cutoff,Editing_cutoff,ASE_cutoff,SNP_cutoff) ## this sequential is important! It determines the label seqeuntial mannually set at melt step!
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
    colnames(Outlier) <- c("Alternative promoter outliers","Expression outliers","Splice outliers","Alternative polyadenyltion outliers","Chimeric RNA","RNA editing outliers","ASE outliers","SNP","UnionN")
    Outlier[is.na(Outlier)] <- 0
    Outlier$Gene <- rownames(Outlier)
    Outlier_plot <- melt(Outlier)
    Outlier_plot$Gene_symbol <- as.character(lapply(strsplit(as.character(Outlier_plot$Gene)," |",fixed = TRUE),function(x) x[1]))
    colnames(Outlier_plot) <- c("Gene","Feature","OutlierN","Gene_symbol")
    
    
    #plot
    #Union
    RNA_Union <- Outlier_plot[grep("UnionN",Outlier_plot$Feature),]
    RNA_Union$Gene <- factor(RNA_Union$Gene,levels=RNA_Union[order(RNA_Union$OutlierN,decreasing = FALSE),]$Gene)
    RNA_Union$Gene_symbol <- factor(RNA_Union$Gene_symbol,levels=RNA_Union[order(RNA_Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
    p_Union_RNA <- ggplot(RNA_Union,aes(x=Gene_symbol,y=OutlierN/sampleN_RNA,fill = Feature))+
      scale_fill_manual(values = c("#138F6A"))+
      geom_bar(stat = "identity")+
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
    RNA_Features$Gene <- factor(RNA_Features$Gene,levels=RNA_Union[order(RNA_Union$OutlierN,decreasing = FALSE),]$Gene)
    RNA_Features$Gene_symbol <- factor(RNA_Features$Gene_symbol,levels=RNA_Union[order(RNA_Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
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
    DNA_nucleosome_cutoff <- DNA_nucleosome_targetgenes[grep(cutoff,colnames(DNA_nucleosome_targetgenes))]
    DNA_MeDIP_cutoff <- DNA_MeDIP_targetgenes[grep(cutoff,colnames(DNA_MeDIP_targetgenes))]
    
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
    DNA_Union$Gene <- factor(DNA_Union$Gene,levels=RNA_Union[order(RNA_Union$OutlierN,decreasing = FALSE),]$Gene)
    DNA_Union$Gene_symbol <- factor(DNA_Union$Gene_symbol,levels=RNA_Union[order(RNA_Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
    p_Union_DNA_rev <- ggplot(DNA_Union,aes(x=Gene_symbol,y=OutlierN/sampleN_DNA,fill = Feature))+
      scale_fill_manual(values = c("#EEC900"))+
      geom_bar(stat = "identity")+
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
        axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))+
      scale_y_reverse(breaks = c(0,0.25,0.50,0.75),labels = c("0.0","","0.5",""),expand = c(0,0),limits = c(0.75,0))
    }
    
    #different features contribution
    DNA_Features <- Outlier_plot[-grep("UnionN",Outlier_plot$Feature),]
    DNA_Features$Gene <- factor(DNA_Features$Gene,levels=DNA_Union[order(DNA_Union$OutlierN,decreasing = FALSE),]$Gene)
    DNA_Features$Gene_symbol <- factor(DNA_Features$Gene_symbol,levels=DNA_Union[order(DNA_Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
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
    DNA_CNV_cutoff <- DNA_CNV_targetgenes[grep(cutoff,colnames(DNA_CNV_targetgenes))]
    DNA_nucleosome_cutoff <- DNA_nucleosome_targetgenes[grep(cutoff,colnames(DNA_nucleosome_targetgenes))]
    DNA_MeDIP_cutoff <- DNA_MeDIP_targetgenes[grep(cutoff,colnames(DNA_MeDIP_targetgenes))]
    
    Alt_cutoff <- Alt_targetgenes[grep(cutoff,colnames(Alt_targetgenes))]
    Expression_cutoff <- Expression_targetgenes[grep(cutoff,colnames(Expression_targetgenes))]
    Splicing_cutoff <- Splicing_targetgenes[grep(cutoff,colnames(Splicing_targetgenes))]
    
    APA_cutoff <- APA_targetgenes[grep(cutoff,colnames(APA_targetgenes))]
    chimeric_cutoff <- chimeric_targetgenes[grep(cutoff,colnames(chimeric_targetgenes))]
    
    Editing_cutoff <- Editing_targetgenes[grep(cutoff,colnames(Editing_targetgenes))]
    ASE_cutoff <- ASE_targetgenes[grep(cutoff,colnames(ASE_targetgenes))]
    SNP_cutoff <- SNP_targetgenes[grep(cutoff,colnames(SNP_targetgenes))]
    
    features <- cbind(DNA_CNV_cutoff,DNA_nucleosome_cutoff,DNA_MeDIP_cutoff,Alt_cutoff,Expression_cutoff,Splicing_cutoff,APA_cutoff,chimeric_cutoff,Editing_cutoff,ASE_cutoff,SNP_cutoff) ## this sequential is important! It determines the label seqeuntial mannually set at melt step!
    
    #melt
    library(reshape)
    Outlier <- features[,grep("_ratio",colnames(features))]
    colnames(Outlier) <- c("DNA: Copy-number outlier","DNA: Nucleosome occupancy outliers","DNA: Methylation outliers","Alternative promoter outliers","Expression outliers","Splice outliers","Alternative polyadenyltion outliers","Chimeric RNA","RNA editing outliers","ASE outliers","RNA SNP")
    #colnames(Outlier) <- c("DNA_CNV","DNA_nucleosome","DNA_MeDIP","Altpromoter","Expression","Splicing","APA","chimeric","Editing","ASE","SNP")
    Outlier[is.na(Outlier)] <- 0
    Outlier$Gene <- rownames(Outlier)
    Outlier_plot <- melt(Outlier)
    Outlier_plot$Gene_symbol <- as.character(lapply(strsplit(as.character(Outlier_plot$Gene)," |",fixed = TRUE),function(x) x[1]))
    colnames(Outlier_plot) <- c("Gene","Feature","OutlierRatio","Gene_symbol")
    
    
    #plot
    #different features contribution
    Features <- Outlier_plot
    Features$Gene <- factor(Features$Gene,levels=RNA_Union[order(RNA_Union$OutlierN,decreasing = FALSE),]$Gene)
    Features$Gene_symbol <- factor(Features$Gene_symbol,levels=RNA_Union[order(RNA_Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
    Features_color <- c("DNA: Copy-number outlier"="#FFE600",
                        "DNA: Nucleosome occupancy outliers"="#60AFFE",
                        "DNA: Methylation outliers"="#CDC9C9",
                        "Alternative promoter outliers"="#E9C2A6",
                        "Expression outliers"="#A5435C",
                        "Splice outliers"="#C5E3BF",
                        "Alternative polyadenyltion outliers"="#003F87",
                        "Chimeric RNA"="#FF3D0D",
                        "RNA editing outliers"="#324F17",
                        "ASE outliers"="#87CEFF",
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
  
  }
}

#KEGG and GO_CC
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/DE")
  GO <-read.table("./plasma/plasma_down_GO_CC.txt",header = T,sep="\t",quote = "")
  GO$Term <- as.character(lapply(strsplit(as.character(GO$Term),"~",fixed = TRUE),function(x) x[2]))
  
  {
    ggplot(data=GO[1:10,])+
      geom_bar(aes(x=reorder(Term,Count),y=Count, fill=-log10(PValue)),
               stat='identity') +
      coord_flip() +
      scale_fill_gradient(expression(-log["10"](P.value)),low="dark blue", high = "red")+
      theme_light()+
      theme(panel.border = element_rect(linetype="solid",,size=0.7,colour = "black"),
            #axis.title.x=element_text(face="bold",size=12),
            #axis.text.x=element_text(face="bold",size=12),
            #axis.text.y=element_text(face="bold",size=12),
            axis.ticks.length = unit(0, "cm"),
            #axis.ticks=element_blank(),
            #axis.line = element_line(colour="black",size=1.2,linetype="solid"),
            #legend.position = "bottom",
            legend.text = element_text(face="bold",size=14),
            legend.direction = NULL,
            legend.title = element_text(face="bold",size=14),
            #legend.key = element_rect(),
            #panel.grid = element_blank(),
            #strip.background = element_rect(fill="white",color="transparent"),
            title = element_text(face="bold",size=24),
            text = element_text(face="bold",size=20),
            axis.text = element_text(face="bold",size=20,colour = "black"))+
      xlab("") +
      ylab("Gene count")
  }
  
  ##KEGG
  KEGG <- read.table("./tissue/tissue_up_KEGG.txt",header=T,sep="\t",quote = "")
  KEGG$Term <- as.character(lapply(strsplit(as.character(KEGG$Term),":",fixed = TRUE),function(x) x[2]))
  ggplot(KEGG[1:10,],aes(x=Fold.Enrichment,y=reorder(Term,Fold.Enrichment)))+
    geom_point(aes(size=Count,color=-1*log10(PValue)))+
    scale_colour_gradient(low="green",high="red")+
    labs(color=expression(-log[10](P.value)),
         size="Gene number",
         x="Fold enrichment")+
    theme_bw()+
    theme(axis.text = element_text(size = rel(2),face="bold",colour = "black"),
          axis.title.x = element_text(size=rel(2),face="bold"),
          axis.title.y = element_blank(),
          legend.text = element_text(size=14),
          legend.title = element_text(face="bold",size=14))
}


#RNA_biotype comparison between specimen
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/02.Quality control")
  specimen_biotype <- read.csv("specimens_RNA_biotype_forplot.csv",header = TRUE, row.names = 1)
  
  specimen_biotype$sncRNA <- specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$srpRNA+specimen_biotype$tRNA+specimen_biotype$Y_RNA
  specimen_biotype$others <- specimen_biotype$tucpRNA+specimen_biotype$unannotated
  
  specimen_biotype_modified <- specimen_biotype[,c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","others")]
  
  specimen_biotype_forplot <- melt(specimen_biotype_modified)
  specimen_biotype_forplot$Type <- rep(rownames(specimen_biotype_modified),ncol(specimen_biotype_modified))
  colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Type")
  
  specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                             levels = c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","others"))
  specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("Plasma-Healthy","Plasma-CRC",
                                                                                 "PBMC-Healthy","PBMC-CRC",
                                                                                 "Tissue-Healthy","Tissue-CRC"))
  ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
    geom_bar(stat = "identity",position = "fill")+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(10,10,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      #axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
      axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
      axis.title.x = element_text(face="bold", color="black", size=20),
      axis.title.y = element_text(face="bold",color="black", size=20))+
    scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01),limits = c(0,1))+
    scale_fill_nejm(alpha = 1)
  
  specimen_sncRNA <- specimen_biotype[,c("snRNA","snoRNA","srpRNA","tRNA","Y_RNA")]
  specimen_sncRNA_forplot <- melt(specimen_sncRNA)
  specimen_sncRNA_forplot$Type <- rep(rownames(specimen_sncRNA),ncol(specimen_sncRNA))
  colnames(specimen_sncRNA_forplot) <- c("Biotype","Ratio","Type")
  specimen_sncRNA_forplot$Biotype <- factor(specimen_sncRNA_forplot$Biotype,
                                             levels = c("srpRNA","Y_RNA","tRNA","snRNA","snoRNA"))
  specimen_sncRNA_forplot$Type <- factor(specimen_sncRNA_forplot$Type,levels=c("Plasma-Healthy","Plasma-CRC",
                                                                                 "PBMC-Healthy","PBMC-CRC",
                                                                                 "Tissue-Healthy","Tissue-CRC"))
  ggplot(specimen_sncRNA_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
    geom_bar(stat = "identity",position = "fill")+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(10,10,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      #axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
      axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
      axis.title.x = element_text(face="bold", color="black", size=20),
      axis.title.y = element_text(face="bold",color="black", size=20))+
    scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01),limits = c(0,1))+
    scale_fill_manual(values = c("#5D478B","#8968CD","#AB82FF","#CCCCFF","#D9D9F3"))
}

#RNA_biotype based on gencode v1
{
#CRC
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/10.RNA_biotypes/20210815_version1/")
specimen_biotype <- read.csv("Specimen_biotypes.csv",header = TRUE, row.names = 1)

specimen_biotype <- specimen_biotype[c("CRC-PKU-27-PBMC","CRC-PKU-28-PBMC","CRC-PKU-29-PBMC","CRC-PKU-30-PBMC","CRC-PKU-32-PBMC",
                                       "CRC-PKU-27-pico","CRC-PKU-28-pico","CRC-PKU-29-pico","CRC-PKU-30-pico","CRC-PKU-32-pico",
                                       "CRC-PKU-27-T","CRC-PKU-28-T","CRC-PKU-29-T","CRC-PKU-30-T","CRC-PKU-32-T"),]

specimen_biotype$MT <- specimen_biotype$MT_tRNA+specimen_biotype$MT_mRNA+specimen_biotype$MT_rRNA+specimen_biotype$MT_lncRNA+specimen_biotype$MT_exon
specimen_biotype$HG <- specimen_biotype$circRNA+specimen_biotype$mRNA+specimen_biotype$lncRNA+specimen_biotype$pseudogene+specimen_biotype$tRNA+specimen_biotype$srpRNA+
                       specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA+specimen_biotype$tucpRNA+specimen_biotype$exon
specimen_biotype$sncRNA <- specimen_biotype$tRNA+specimen_biotype$srpRNA+specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA+
                           specimen_biotype$tucpRNA
specimen_biotype$non_exon <- specimen_biotype$repeats+specimen_biotype$promoter+specimen_biotype$enhancer

specimen_biotype_all <- specimen_biotype[,c("MT","HG")]
specimen_biotype_MT <- specimen_biotype[,c("MT_tRNA","MT_mRNA","MT_rRNA","MT_lncRNA")]
specimen_biotype_HG <- specimen_biotype[,c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","tucpRNA")]
specimen_biotype_sncRNA <- specimen_biotype[,c("tRNA","srpRNA","snoRNA","snRNA","Y_RNA")]
specimen_biotype_yunfan <- specimen_biotype[,c("circRNA","lncRNA","MT_rRNA","MT_tRNA","MT_mRNA","mRNA","exon","non_exon")]

specimen_biotype_all_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_all),margin = 1))
specimen_biotype_MT_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_MT),margin = 1))
specimen_biotype_HG_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_HG),margin = 1))
specimen_biotype_sncRNA_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_sncRNA),margin = 1))
specimen_biotype_yunfan_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_yunfan),margin = 1))

#yunfan
{
  specimen_biotype_forplot <- melt(specimen_biotype_yunfan_ratio)
  specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_yunfan_ratio),ncol(specimen_biotype_yunfan_ratio))
  specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_yunfan_ratio))
  colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
  
  specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                             levels = c("circRNA","lncRNA","MT_rRNA","MT_tRNA","MT_mRNA","mRNA","exon","non_exon"))
  specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                 "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                 "Tissue-Healthy","Tissue-CRC"))
  #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
  yunfan <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
    geom_bar(stat = "identity",position = "fill")+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(
      #legend.position="none",
      plot.margin = unit(x=c(10,10,10,10),units="pt"),
      legend.position="top",
      legend.direction="vertical",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      legend.title = element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      #axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
      axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
      axis.title.x = element_text(face="bold", color="black", size=20),
      axis.title.y = element_text(face="bold",color="black", size=20))+
    scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
    scale_fill_nejm()
}

#MT+HG
{
specimen_biotype_forplot <- melt(specimen_biotype_all)
specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_all),ncol(specimen_biotype_all))
specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_all))
colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")

specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                           levels = c("MT","repeats","promoter","enhancer","others", "HG"))
specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                               "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                               "Tissue-Healthy","Tissue-CRC"))
#specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
all <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
  geom_bar(stat = "identity",position = "fill")+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(
    #legend.position="none",
    plot.margin = unit(x=c(10,10,10,10),units="pt"),
    legend.position="top",
    legend.direction="vertical",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    legend.title = element_blank(),
    #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    #axis.ticks = element_blank(),
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
    axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
    axis.title.x = element_text(face="bold", color="black", size=20),
    axis.title.y = element_text(face="bold",color="black", size=20))+
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
  scale_fill_manual(values = c("#BE2625","#A0522D","#EE9A00","#006633","#8C8C8C","#5D478B"))
}
#MT
{
  specimen_biotype_forplot <- melt(specimen_biotype_MT)
  specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_MT),ncol(specimen_biotype_MT))
  specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_MT))
  colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
  
  specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                             levels = c("MT_tRNA","MT_rRNA","MT_lncRNA","MT_exon","MT_mRNA"))
  specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                 "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                 "Tissue-Healthy","Tissue-CRC"))
  #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
  MT <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
    geom_bar(stat = "identity",position = "fill")+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(
      #legend.position="none",
      plot.margin = unit(x=c(10,10,10,10),units="pt"),
      legend.position="top",
      legend.direction="vertical",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      legend.title = element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      #axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
      axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
      axis.title.x = element_text(face="bold", color="black", size=20),
      axis.title.y = element_text(face="bold",color="black", size=20))+
    scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
    scale_fill_manual(values = c("#B22222","#CC3232","#EE3B3B","#FF6666","#FFCCCC"))
}
#HG
{
  specimen_biotype_forplot <- melt(specimen_biotype_HG)
  specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_HG),ncol(specimen_biotype_HG))
  specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_HG))
  colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
  
  specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                             levels = c("circRNA","lncRNA","pseudogene","sncRNA","tucpRNA","mRNA"))
  specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                 "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                 "Tissue-Healthy","Tissue-CRC"))
  specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
  HG <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
    geom_bar(stat = "identity",position = "fill")+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(
      #legend.position="none",
      plot.margin = unit(x=c(10,10,10,10),units="pt"),
      legend.position="top",
      legend.direction="vertical",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      legend.title = element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      #axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
      axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
      axis.title.x = element_text(face="bold", color="black", size=20),
      axis.title.y = element_text(face="bold",color="black", size=20))+
    scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
    scale_fill_manual(values = c("#5D478B","#8968CD","#AB82FF","#CCCCFF","#D9D9F3","#F8F8FF"))
}
#sncRNA
{
  specimen_biotype_forplot <- melt(specimen_biotype_sncRNA)
  specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_sncRNA),ncol(specimen_biotype_sncRNA))
  specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_sncRNA))
  colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
  
  specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                             levels = c("Y_RNA","tRNA","srpRNA","snoRNA","snRNA"))
  specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                 "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                 "Tissue-Healthy","Tissue-CRC"))
  #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
  sncRNA <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
    geom_bar(stat = "identity",position = "fill")+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(
      #legend.position="none",
      plot.margin = unit(x=c(10,10,10,10),units="pt"),
      legend.position="top",
      legend.direction="vertical",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      legend.title = element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      #axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
      axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
      axis.title.x = element_text(face="bold", color="black", size=20),
      axis.title.y = element_text(face="bold",color="black", size=20))+
    scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
    scale_fill_manual(values = c("#00688B","#65909A","#50A6C2","#009ACD","#ADD8E6","#F0F8FF"))
}
ggarrange(all,HG,sncRNA,align = "h", ncol = 3, nrow =1)
}
#NC
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/10.RNA_biotypes/20210815_version1/")
  specimen_biotype <- read.csv("Specimen_biotypes.csv",header = TRUE, row.names = 1)
  
  specimen_biotype <- specimen_biotype[c("NC-PKU-9-PBMC","NC-PKU-10-PBMC","NC-PKU-11-PBMC","NC-PKU-12-PBMC","NC-PKU-13-PBMC","NC-PKU-14-PBMC",
                                         "NC-PKU-9-pico","NC-PKU-10-pico","NC-PKU-11-pico","NC-PKU-12-pico","NC-PKU-13-pico","NC-PKU-14-pico",
                                         "CRC-PKU-27-N","CRC-PKU-28-N","CRC-PKU-29-N","CRC-PKU-30-N","CRC-PKU-32-N"),]
  
  specimen_biotype$MT <- specimen_biotype$MT_tRNA+specimen_biotype$MT_mRNA+specimen_biotype$MT_rRNA+specimen_biotype$MT_lncRNA+specimen_biotype$MT_exon
  specimen_biotype$HG <- specimen_biotype$circRNA+specimen_biotype$mRNA+specimen_biotype$lncRNA+specimen_biotype$pseudogene+specimen_biotype$tRNA+specimen_biotype$srpRNA+
    specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA+specimen_biotype$tucpRNA+specimen_biotype$exon
  specimen_biotype$sncRNA <- specimen_biotype$tRNA+specimen_biotype$srpRNA+specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA+
    specimen_biotype$tucpRNA
  
  specimen_biotype_all <- specimen_biotype[,c("MT","HG","repeats","promoter","enhancer","others")]
  specimen_biotype_MT <- specimen_biotype[,c("MT_tRNA","MT_mRNA","MT_rRNA","MT_lncRNA","MT_exon")]
  specimen_biotype_HG <- specimen_biotype[,c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","exon")]
  specimen_biotype_sncRNA <- specimen_biotype[,c("tRNA","srpRNA","snoRNA","snRNA","Y_RNA","tucpRNA")]
  
  #MT+HG
  {
    specimen_biotype_forplot <- melt(specimen_biotype_all)
    specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_all),ncol(specimen_biotype_all))
    specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_all))
    colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
    
    specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                               levels = c("MT","repeats","promoter","enhancer","others", "HG"))
    specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                   "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                   "Tissue-Healthy","Tissue-CRC"))
    #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
    all <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
      geom_bar(stat = "identity",position = "fill")+
      xlab("")+
      ylab("")+
      theme_bw()+
      theme(
        #legend.position="none",
        plot.margin = unit(x=c(10,10,10,10),units="pt"),
        legend.position="top",
        legend.direction="vertical",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.title = element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        #axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))+
      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
      scale_fill_manual(values = c("#BE2625","#A0522D","#EE9A00","#006633","#8C8C8C","#5D478B"))
  }
  
  #MT
  {
    specimen_biotype_forplot <- melt(specimen_biotype_MT)
    specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_MT),ncol(specimen_biotype_MT))
    specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_MT))
    colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
    
    specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                               levels = c("MT_tRNA","MT_rRNA","MT_lncRNA","MT_exon","MT_mRNA"))
    specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                   "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                   "Tissue-Healthy","Tissue-CRC"))
    #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
    MT <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
      geom_bar(stat = "identity",position = "fill")+
      xlab("")+
      ylab("")+
      theme_bw()+
      theme(
        #legend.position="none",
        plot.margin = unit(x=c(10,10,10,10),units="pt"),
        legend.position="top",
        legend.direction="vertical",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.title = element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        #axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))+
      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
      scale_fill_manual(values = c("#B22222","#CC3232","#EE3B3B","#FF6666","#FFCCCC"))
  }
  
  #HG
  {
    specimen_biotype_forplot <- melt(specimen_biotype_HG)
    specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_HG),ncol(specimen_biotype_HG))
    specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_HG))
    colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
    
    specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                               levels = c("circRNA","lncRNA","pseudogene","sncRNA","exon","mRNA"))
    specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                   "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                   "Tissue-Healthy","Tissue-CRC"))
    specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
    HG <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
      geom_bar(stat = "identity",position = "fill")+
      xlab("")+
      ylab("")+
      theme_bw()+
      theme(
        #legend.position="none",
        plot.margin = unit(x=c(10,10,10,10),units="pt"),
        legend.position="top",
        legend.direction="vertical",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.title = element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        #axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))+
      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
      scale_fill_manual(values = c("#5D478B","#8968CD","#AB82FF","#CCCCFF","#D9D9F3","#F8F8FF"))
  }
  #sncRNA
  {
    specimen_biotype_forplot <- melt(specimen_biotype_sncRNA)
    specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_sncRNA),ncol(specimen_biotype_sncRNA))
    specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_sncRNA))
    colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
    
    specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                               levels = c("Y_RNA","tRNA","srpRNA","snoRNA","snRNA","tucpRNA"))
    specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                   "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                   "Tissue-Healthy","Tissue-CRC"))
    #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
    sncRNA <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
      geom_bar(stat = "identity",position = "fill")+
      xlab("")+
      ylab("")+
      theme_bw()+
      theme(
        #legend.position="none",
        plot.margin = unit(x=c(10,10,10,10),units="pt"),
        legend.position="top",
        legend.direction="vertical",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.title = element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        #axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))+
      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
      scale_fill_manual(values = c("#00688B","#65909A","#50A6C2","#009ACD","#ADD8E6","#F0F8FF"))
  }
  ggarrange(MT,all,HG,sncRNA,align = "h", ncol = 4, nrow =1)
}
#comparison between specimen in paired CRC
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/10.RNA_biotypes/20210815_version1/")
  specimen_biotype <- read.csv("Specimen_biotypes.csv",header = TRUE, row.names = 1)
  
  specimen_biotype <- specimen_biotype[c("CRC-PKU-27-PBMC","CRC-PKU-28-PBMC","CRC-PKU-29-PBMC","CRC-PKU-30-PBMC","CRC-PKU-32-PBMC",
                                         "CRC-PKU-27-pico","CRC-PKU-28-pico","CRC-PKU-29-pico","CRC-PKU-30-pico","CRC-PKU-32-pico",
                                         "CRC-PKU-27-T","CRC-PKU-28-T","CRC-PKU-29-T","CRC-PKU-30-T","CRC-PKU-32-T"),]
  
  specimen_biotype$MT <- specimen_biotype$MT_tRNA+specimen_biotype$MT_mRNA+specimen_biotype$MT_rRNA+specimen_biotype$MT_lncRNA+specimen_biotype$MT_exon
  specimen_biotype$HG <- specimen_biotype$circRNA+specimen_biotype$mRNA+specimen_biotype$lncRNA+specimen_biotype$pseudogene+specimen_biotype$tRNA+specimen_biotype$srpRNA+
    specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA+specimen_biotype$tucpRNA+specimen_biotype$exon
  specimen_biotype$sncRNA <- specimen_biotype$tRNA+specimen_biotype$srpRNA+specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA+
    specimen_biotype$tucpRNA
  
  specimen_biotype_all <- specimen_biotype[,c("MT","HG")]
  specimen_biotype_MT <- specimen_biotype[,c("MT_tRNA","MT_mRNA","MT_rRNA","MT_lncRNA","MT_exon")]
  specimen_biotype_HG <- specimen_biotype[,c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","tucpRNA")]
  specimen_biotype_sncRNA <- specimen_biotype[,c("tRNA","srpRNA","snoRNA","snRNA","Y_RNA")]
  
  specimen_biotype_all_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_all),margin = 1))
  specimen_biotype_MT_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_MT),margin = 1))
  specimen_biotype_HG_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_HG),margin = 1))
  specimen_biotype_sncRNA_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_sncRNA),margin = 1))
  
  #MT
  specimen_biotype_forplot <- melt(specimen_biotype_MT_ratio)
  specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_MT_ratio),ncol(specimen_biotype_MT_ratio))
  specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_MT_ratio))
  specimen_biotype_forplot$ID <- as.character(lapply(strsplit(specimen_biotype_forplot$Sample,"-"),function(x) paste(x[1],x[2],x[3], sep ="-", collapse =NULL)))
  colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type","ID")
  
  #all
  #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
  #                                           levels = c("MT", "HG"))
  #HG
  #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
  #                                           levels = c("circRNA","lncRNA","pseudogene","sncRNA","tucpRNA","mRNA"))
  #sncRNA
  #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
  #                                           levels = c("Y_RNA","tRNA","srpRNA","snoRNA","snRNA"))
  #MT
  specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                             levels = c("MT_tRNA","MT_rRNA","MT_lncRNA","MT_exon","MT_mRNA"))
  
  specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                 "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                 "Tissue-Healthy","Tissue-CRC"))
  my_comparisons <- list(c("Plasma-CRC","PBMC-CRC"),c("Plasma-CRC","Tissue-CRC"))
  p <- ggplot(specimen_biotype_forplot[which(specimen_biotype_forplot$Biotype=="MT_tRNA"),],aes(x=Type,y=Ratio,fill=Type))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
    geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
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
  
  m=10
  if(m>0){
    p <- p+geom_line(aes(group = ID))+
      stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test",
                              method.args = list(alternative = "greater",paired = TRUE),
                              label = "p.signif"
    )+labs(x="",y="Rario (%)",title=paste0("MT_tRNA.wilcox.greater"), face="bold",fill="Specimen")
  } else {
    p <- p+geom_line(aes(group = ID))+
      stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test",
                              method.args = list(alternative = "less",paired = TRUE),
                              label = "p.signif"
    )+labs(x="",y="Rario (%)",title=paste0("MT_tRNA.wilcox.test.less"), face="bold",fill="Batch")
  }
}

#comparison between CRC and NC in patients and healthy donor
{
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/10.RNA_biotypes/20210815_version1/")
    specimen_biotype <- read.csv("Specimen_biotypes.csv",header = TRUE, row.names = 1)
    
    specimen_biotype <- specimen_biotype[c("CRC-PKU-27-PBMC","CRC-PKU-28-PBMC","CRC-PKU-29-PBMC","CRC-PKU-30-PBMC","CRC-PKU-31-PBMC","CRC-PKU-32-PBMC",
                                           "NC-PKU-9-PBMC","NC-PKU-10-PBMC","NC-PKU-11-PBMC","NC-PKU-12-PBMC","NC-PKU-13-PBMC","NC-PKU-14-PBMC",
                                           "CRC-2280766",
                                           "CRC-2327118",
                                           "CRC-2334660",
                                           "CRC-2362297",
                                           "CRC-2372270",
                                           "CRC-2380493",
                                           "CRC-2382513",
                                           "CRC-2383018",
                                           "CRC-2383387",
                                           "CRC-2383420",
                                           "CRC-2383555",
                                           "CRC-2383730",
                                           "CRC-2384058",
                                           "CRC-2384076",
                                           "CRC-2384108",
                                           "CRC-2384226",
                                           "CRC-2384272",
                                           "CRC-2384377",
                                           "CRC-2384422",
                                           "CRC-2384541",
                                           "CRC-2384548",
                                           "CRC-2384796",
                                           "CRC-2386412",
                                           "CRC-2386612",
                                           "CRC-2387327",
                                           "CRC-2388859",
                                           "CRC-2388932",
                                           "CRC-2389212",
                                           "CRC-2390502",
                                           "CRC-2390539",
                                           "CRC-2390716",
                                           "CRC-2391097",
                                           "CRC-2393807",
                                           "CRC-2396217",
                                           "CRC-2396837",
                                           "CRC-2397180",
                                           "CRC-2398897",
                                           "CRC-2399129",
                                           "CRC-2405108",
                                           "CRC-2406345",
                                           "CRC-2406552",
                                           "CRC-2407821",
                                           "CRC-2408137",
                                           "CRC-2408162",
                                           "CRC-2408206",
                                           "CRC-2408346",
                                           "CRC-2408784",
                                           "CRC-2410759",
                                           "CRC-2410966",
                                           "CRC-2411883",
                                           "CRC-2411896",
                                           "CRC-2412048",
                                           "CRC-2412066",
                                           "CRC-9299429",
                                           "CRC-PKU-10-pico",
                                           "CRC-PKU-12-pico",
                                           "CRC-PKU-13-pico",
                                           "CRC-PKU-14-pico",
                                           "CRC-PKU-15-pico",
                                           "CRC-PKU-16-pico",
                                           "CRC-PKU-17-pico",
                                           "CRC-PKU-18-pico",
                                           "CRC-PKU-19-pico",
                                           "CRC-PKU-20-pico",
                                           "CRC-PKU-21-pico",
                                           "CRC-PKU-22-pico",
                                           "CRC-PKU-23-pico",
                                           "CRC-PKU-24-pico",
                                           "CRC-PKU-25-pico",
                                           "CRC-PKU-26-pico",
                                           "CRC-PKU-27-pico",
                                           "CRC-PKU-28-pico",
                                           "CRC-PKU-29-pico",
                                           "CRC-PKU-2-pico",
                                           "CRC-PKU-30-pico",
                                           "CRC-PKU-31-pico",
                                           "CRC-PKU-32-pico",
                                           "CRC-PKU-33-pico",
                                           "CRC-PKU-34-pico",
                                           "CRC-PKU-35-pico",
                                           "CRC-PKU-36-pico",
                                           "CRC-PKU-37-pico",
                                           "CRC-PKU-38-pico",
                                           "CRC-PKU-39-pico",
                                           "CRC-PKU-3-pico",
                                           "CRC-PKU-40-pico",
                                           "CRC-PKU-41-pico",
                                           "CRC-PKU-4-pico",
                                           "CRC-PKU-5-pico",
                                           "CRC-PKU-6-pico",
                                           "CRC-PKU-8-pico",
                                           "CRC-PKU-9-pico",
                                           "CRC-PKU-mix1-1-pico",
                                           "CRC-PKU-mix1-2-pico",
                                           "CRC-PKU-mix1-3-pico",
                                           "NC-PKU-10-pico",
                                           "NC-PKU-11-pico",
                                           "NC-PKU-12-pico",
                                           "NC-PKU-13-pico",
                                           "NC-PKU-14-pico",
                                           "NC-PKU-2359973",
                                           "NC-PKU-2392860",
                                           "NC-PKU-2392931",
                                           "NC-PKU-2392932",
                                           "NC-PKU-2392993",
                                           "NC-PKU-2394780",
                                           "NC-PKU-2396996",
                                           "NC-PKU-2397512",
                                           "NC-PKU-2397878",
                                           "NC-PKU-2399728",
                                           "NC-PKU-28",
                                           "NC-PKU-33",
                                           "NC-PKU-34",
                                           "NC-PKU-37",
                                           "NC-PKU-48",
                                           "NC-PKU-50",
                                           "NC-PKU-56",
                                           "NC-PKU-69",
                                           "NC-PKU-8111237",
                                           "NC-PKU-9-pico",
                                           "NC-PKU-mix1-1-pico",
                                           "NC-PKU-mix1-2-pico",
                                           "NC-PKU-mix1-3-pico",
                                           "NC-PKU-mix15-pico",
                                           "NC-PKU-mix16-pico",
                                           "NC-PKU-mix17-pico",
                                           "NC-PKU-mix18-pico",
                                           "NC-PKU-mix19-pico",
                                           
                                           "NC-PKU-mix20-pico",
                                           "NC-PKU-mix2-1-pico",
                                           "NC-PKU-mix21-pico",
                                           "NC-PKU-mix2-2-pico",
                                           "NC-PKU-mix22-pico",
                                           "NC-PKU-mix24-pico",
                                           "NC-PKU-mix25-pico",
                                           "NC-PKU-mix26-pico",
                                           "NC-PKU-mix27-pico",
                                           "NC-PKU-mix28-pico",
                                           "NC-PKU-mix29-pico",
                                           
                                           "NC-PKU-mix30-pico",
                                           "NC-PKU-mix3-1-pico",
                                           "NC-PKU-mix3-2-pico",
                                           "NC-PKU-mix3-3-pico",
                                           "NC-PKU-mix3-4-pico",
                                           "NC-PKU-mix3-5-pico",
                                         
                                           "NC-PKU-mix4-1-pico",
                                           "NC-PKU-mix4-2-pico",
                                           
                                           "NC-PKU-mix5-1-pico",
                                           "NC-PKU-mix5-2-pico",
                                           "NC-PKU-mix5-3-pico",
                                           
                                           "NC-PKU-mix6-1-pico",
                                           "NC-PKU-mix6-2-pico",
                                           "NC-PKU-mix6-3-pico",
                                           "NC-PKU-mix6-4-pico",
                                           "NC-PKU-mix6-5-pico",
                                           "NC-PKU-mix6-6-pico",
                                           
                                           "NC-PKU-mix7-1-pico",
                                           "NC-PKU-mix7-2-pico",
                                           
                                           "NC-PKU-mix8-1-pico",
                                           "NC-PKU-mix8-2-pico",
                                           "NC-PKU-mix8-3-pico",
                                           "NC-PKU-mix8-4-pico",
                                           
                                           "NC-PKU-mix9-1-pico",
                                           "NC-PKU-mix9-2-pico",
                                           "NC-PKU-mix9-pico",
                                           
                                           #"CRC-PKU-27-pico","CRC-PKU-28-pico","CRC-PKU-29-pico","CRC-PKU-30-pico","CRC-PKU-31-pico","CRC-PKU-32-pico","CRC-PKU-34-pico","CRC-PKU-35-pico","CRC-PKU-36-pico","CRC-PKU-37-pico","CRC-PKU-38-pico","CRC-PKU-39-pico","CRC-PKU-40-pico","CRC-PKU-41-pico",
                                           #"NC-PKU-9-pico","NC-PKU-10-pico","NC-PKU-11-pico","NC-PKU-12-pico","NC-PKU-13-pico","NC-PKU-14-pico",
                                           "CRC-PKU-27-T","CRC-PKU-28-T","CRC-PKU-29-T","CRC-PKU-30-T","CRC-PKU-32-T","CRC-PKU-34-T","CRC-PKU-35-T","CRC-PKU-36-T","CRC-PKU-37-T","CRC-PKU-38-T","CRC-PKU-39-T","CRC-PKU-40-T","CRC-PKU-41-T",
                                           "CRC-PKU-27-N","CRC-PKU-28-N","CRC-PKU-29-N","CRC-PKU-30-N","CRC-PKU-32-N","CRC-PKU-34-N","CRC-PKU-35-N","CRC-PKU-36-N","CRC-PKU-37-N","CRC-PKU-38-N","CRC-PKU-39-N","CRC-PKU-40-N","CRC-PKU-41-N"),]
    
    specimen_biotype$MT <- specimen_biotype$MT_tRNA+specimen_biotype$MT_mRNA+specimen_biotype$MT_rRNA+specimen_biotype$MT_lncRNA+specimen_biotype$MT_exon
    specimen_biotype$HG <- specimen_biotype$circRNA+specimen_biotype$mRNA+specimen_biotype$lncRNA+specimen_biotype$pseudogene+specimen_biotype$tRNA+specimen_biotype$srpRNA+
      specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA+specimen_biotype$tucpRNA+specimen_biotype$exon
    specimen_biotype$sncRNA <- specimen_biotype$tRNA+specimen_biotype$srpRNA+specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA+
      specimen_biotype$tucpRNA
    
    specimen_biotype_all <- specimen_biotype[,c("MT","HG","repeats","promoter","enhancer","others")]
    specimen_biotype_MT <- specimen_biotype[,c("MT_tRNA","MT_mRNA","MT_rRNA","MT_lncRNA","MT_exon")]
    specimen_biotype_HG <- specimen_biotype[,c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","exon")]
    specimen_biotype_sncRNA <- specimen_biotype[,c("tRNA","srpRNA","snoRNA","snRNA","Y_RNA","tucpRNA")]
    
    specimen_biotype_all_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_all),margin = 1))
    specimen_biotype_MT_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_MT),margin = 1))
    specimen_biotype_HG_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_HG),margin = 1))
    specimen_biotype_sncRNA_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_sncRNA),margin = 1))
    
    #sncRNA
    specimen_biotype_forplot <- melt(specimen_biotype_sncRNA_ratio)
    specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_sncRNA_ratio),ncol(specimen_biotype_sncRNA_ratio))
    specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_sncRNA_ratio))
    specimen_biotype_forplot$ID <- as.character(lapply(strsplit(specimen_biotype_forplot$Sample,"-"),function(x) paste(x[1],x[2],x[3], sep ="-", collapse =NULL)))
    specimen_biotype_forplot$ID <- paste(specimen_biotype_forplot$ID,specimen_biotype_forplot$Type,sep ="-", collapse =NULL)
    specimen_biotype_forplot$ID <- as.character(lapply(strsplit(specimen_biotype_forplot$ID,"-"),function(x) paste(x[1],x[2],x[3],x[4], sep ="-", collapse =NULL)))
    colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type","ID")
    
    #all
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("MT","repeats","promoter","enhancer","others", "HG"))
    #HG
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("circRNA","lncRNA","pseudogene","sncRNA","exon","mRNA"))
    #sncRNA
    specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                               levels = c("Y_RNA","tRNA","srpRNA","snoRNA","snRNA","tucpRNA"))
    #MT
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("MT_tRNA","MT_rRNA","MT_lncRNA","MT_exon","MT_mRNA"))
    
    specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-CRC","PBMC-Healthy",
                                                                                   "Plasma-CRC","Plasma-Healthy","Plasma-STAD",
                                                                                   "Tissue-CRC","Tissue-Healthy"))
    my_comparisons <- list(c("PBMC-CRC","PBMC-Healthy"),c("Plasma-CRC","Plasma-Healthy"),c("Tissue-CRC","Tissue-Healthy"))
    p <- ggplot(specimen_biotype_forplot[which(specimen_biotype_forplot$Biotype=="srpRNA"),],aes(x=Type,y=Ratio,fill=Type))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("red","blue","red","blue","red","blue")) +
      geom_vline(xintercept = c(2.5,4.5),linetype = "dashed",color="grey")+
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
    if(m==0) {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "two.sided",paired = TRUE,exact = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Rario (%)",title=paste0("srpRNA.wilcox"), face="bold",fill="Specimen")
    } else if(m>0) {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "greater",paired = TRUE,exact = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Rario (%)",title=paste0("srpRNA.wilcox.greater"), face="bold",fill="Specimen")
    } else {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "less",paired = TRUE,exact = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Rario (%)",title=paste0("srpRNA.wilcox.test.less"), face="bold",fill="Specimen")
    }
}
  
}

#RNA_biotype based on gencode v2
{
  #CRC
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/10.RNA_biotypes/20210820_version2/")
    specimen_biotype <- read.csv("Specimen_biotypes.csv",header = TRUE, row.names = 1)
    
    specimen_biotype <- specimen_biotype[c("CRC-PKU-27-PBMC","CRC-PKU-28-PBMC","CRC-PKU-29-PBMC","CRC-PKU-30-PBMC","CRC-PKU-32-PBMC",
                                           "CRC-PKU-27-pico","CRC-PKU-28-pico","CRC-PKU-29-pico","CRC-PKU-30-pico","CRC-PKU-32-pico",
                                           "CRC-PKU-27-T","CRC-PKU-28-T","CRC-PKU-29-T","CRC-PKU-30-T","CRC-PKU-32-T"),]
    
    specimen_biotype$MT <- specimen_biotype$MT_tRNA+specimen_biotype$MT_mRNA+specimen_biotype$MT_rRNA+specimen_biotype$MT_lncRNA+specimen_biotype$MT_exon
    specimen_biotype$sncRNA <- specimen_biotype$tRNA+specimen_biotype$srpRNA+specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA
    specimen_biotype$cleaned <- specimen_biotype$genome+specimen_biotype$circRNA+specimen_biotype$unmapped
    
    specimen_biotype$HG <- specimen_biotype$circRNA+specimen_biotype$genome-specimen_biotype$MT
    specimen_biotype$unassigned <- specimen_biotype$unmapped-specimen_biotype$microbe
    
    specimen_biotype_all <- specimen_biotype[,c("MT","HG","microbe","unassigned")]
    specimen_biotype_MT <- specimen_biotype[,c("MT_tRNA","MT_mRNA","MT_rRNA","MT_lncRNA")]
    specimen_biotype_HG <- specimen_biotype[,c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","tucpRNA","misc")]
    specimen_biotype_sncRNA <- specimen_biotype[,c("tRNA","srpRNA","snoRNA","snRNA","Y_RNA")]
    
    specimen_biotype_all_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_all),margin = 1))
    specimen_biotype_MT_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_MT),margin = 1))
    specimen_biotype_HG_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_HG),margin = 1))
    specimen_biotype_sncRNA_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_sncRNA),margin = 1))
    write.csv(specimen_biotype_all_ratio,"all_ratio.csv")
    
    #MT+HG+microbe+unassiganed
    {
      specimen_biotype_forplot <- melt(specimen_biotype_all)
      specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_all),ncol(specimen_biotype_all))
      specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_all))
      colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
      
      specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                                 levels = c("MT","HG","microbe","unassigned"))
      specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                     "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                     "Tissue-Healthy","Tissue-CRC"))
      #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
      all <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
        geom_bar(stat = "identity",position = "fill")+
        xlab("")+
        ylab("")+
        theme_bw()+
        theme(
          #legend.position="none",
          plot.margin = unit(x=c(10,10,10,10),units="pt"),
          legend.position="top",
          legend.direction="vertical",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.title = element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          #axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))+
        scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
        scale_fill_manual(values = c("#BE2625","#A0522D","#EE9A00","#006633","#8C8C8C","#5D478B"))
    }
    #MT
    {
      specimen_biotype_forplot <- melt(specimen_biotype_MT)
      specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_MT),ncol(specimen_biotype_MT))
      specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_MT))
      colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
      
      specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                                 levels = c("MT_tRNA","MT_rRNA","MT_lncRNA","MT_exon","MT_mRNA"))
      specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                     "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                     "Tissue-Healthy","Tissue-CRC"))
      #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
      MT <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
        geom_bar(stat = "identity",position = "fill")+
        xlab("")+
        ylab("")+
        theme_bw()+
        theme(
          #legend.position="none",
          plot.margin = unit(x=c(10,10,10,10),units="pt"),
          legend.position="top",
          legend.direction="vertical",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.title = element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          #axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))+
        scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
        scale_fill_manual(values = c("#B22222","#CC3232","#EE3B3B","#FF6666","#FFCCCC"))
    }
    #HG
    {
      specimen_biotype_forplot <- melt(specimen_biotype_HG)
      specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_HG),ncol(specimen_biotype_HG))
      specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_HG))
      colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
      
      specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                                 levels = c("circRNA","lncRNA","pseudogene","sncRNA","tucpRNA","misc","mRNA"))
      specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                     "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                     "Tissue-Healthy","Tissue-CRC"))
      specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
      HG <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
        geom_bar(stat = "identity",position = "fill")+
        xlab("")+
        ylab("")+
        theme_bw()+
        theme(
          #legend.position="none",
          plot.margin = unit(x=c(10,10,10,10),units="pt"),
          legend.position="top",
          legend.direction="vertical",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.title = element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          #axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))+
        scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
        #scale_fill_manual(values = c("#5D478B","#8968CD","#AB82FF","#CCCCFF","#D9D9F3","#F8F8FF"))
        #scale_fill_aaas()
        scale_fill_manual(values=c("#631779","#ee2200","#008b45","#3b4992","#008280","#bb1921","#5f559b"))

    }
    #sncRNA
    {
      specimen_biotype_forplot <- melt(specimen_biotype_sncRNA)
      specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_sncRNA),ncol(specimen_biotype_sncRNA))
      specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_sncRNA))
      colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
      
      specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                                 levels = c("Y_RNA","tRNA","srpRNA","snoRNA","snRNA"))
      specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                     "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                     "Tissue-Healthy","Tissue-CRC"))
      #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
      sncRNA <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
        geom_bar(stat = "identity",position = "fill")+
        xlab("")+
        ylab("")+
        theme_bw()+
        theme(
          #legend.position="none",
          plot.margin = unit(x=c(10,10,10,10),units="pt"),
          legend.position="top",
          legend.direction="vertical",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.title = element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          #axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))+
        scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
        scale_fill_manual(values = c("#00688B","#65909A","#50A6C2","#009ACD","#ADD8E6","#F0F8FF"))
    }
    ggarrange(MT,all,HG,sncRNA,align = "h", ncol = 4, nrow =1)
  }
  #comparison between specimen in paired CRC
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/10.RNA_biotypes/20210820_version2/")
    specimen_biotype <- read.csv("Specimen_biotypes.csv",header = TRUE, row.names = 1)
    
    specimen_biotype <- specimen_biotype[c("CRC-PKU-27-PBMC","CRC-PKU-28-PBMC","CRC-PKU-29-PBMC","CRC-PKU-30-PBMC","CRC-PKU-32-PBMC",
                                           "CRC-PKU-27-pico","CRC-PKU-28-pico","CRC-PKU-29-pico","CRC-PKU-30-pico","CRC-PKU-32-pico",
                                           "CRC-PKU-27-T","CRC-PKU-28-T","CRC-PKU-29-T","CRC-PKU-30-T","CRC-PKU-32-T"),]
    
    specimen_biotype$MT <- specimen_biotype$MT_tRNA+specimen_biotype$MT_mRNA+specimen_biotype$MT_rRNA+specimen_biotype$MT_lncRNA+specimen_biotype$MT_exon
    specimen_biotype$sncRNA <- specimen_biotype$tRNA+specimen_biotype$srpRNA+specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA
    specimen_biotype$cleaned <- specimen_biotype$genome+specimen_biotype$circRNA+specimen_biotype$unmapped
    
    specimen_biotype$HG <- specimen_biotype$circRNA+specimen_biotype$genome-specimen_biotype$MT
    specimen_biotype$unassigned <- specimen_biotype$unmapped-specimen_biotype$microbe
    
    specimen_biotype_all <- specimen_biotype[,c("MT","HG","microbe","unassigned")]
    specimen_biotype_MT <- specimen_biotype[,c("MT_tRNA","MT_mRNA","MT_rRNA","MT_lncRNA")]
    specimen_biotype_HG <- specimen_biotype[,c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","tucpRNA","misc")]
    specimen_biotype_sncRNA <- specimen_biotype[,c("tRNA","srpRNA","snoRNA","snRNA","Y_RNA")]
    
    specimen_biotype_all_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_all),margin = 1))
    specimen_biotype_MT_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_MT),margin = 1))
    specimen_biotype_HG_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_HG),margin = 1))
    specimen_biotype_sncRNA_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_sncRNA),margin = 1))
    
    #MT
    specimen_biotype_forplot <- melt(specimen_biotype_MT_ratio)
    specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_MT_ratio),ncol(specimen_biotype_MT_ratio))
    specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_MT_ratio))
    specimen_biotype_forplot$ID <- as.character(lapply(strsplit(specimen_biotype_forplot$Sample,"-"),function(x) paste(x[1],x[2],x[3], sep ="-", collapse =NULL)))
    colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type","ID")
    
    #all
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("MT","HG","microbe","unassigned"))
    #HG
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("circRNA","lncRNA","pseudogene","sncRNA","tucpRNA","mRNA"))
    #sncRNA
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("Y_RNA","tRNA","srpRNA","snoRNA","snRNA"))
    #MT
    specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                               levels = c("MT_tRNA","MT_rRNA","MT_lncRNA","MT_exon","MT_mRNA"))
    
    specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                   "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                   "Tissue-Healthy","Tissue-CRC"))
    my_comparisons <- list(c("Plasma-CRC","PBMC-CRC"),c("Plasma-CRC","Tissue-CRC"))
    p <- ggplot(specimen_biotype_forplot[which(specimen_biotype_forplot$Biotype=="MT_rRNA"),],aes(x=Type,y=Ratio,fill=Type))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
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
    
    m=0
    if(m==0) {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "two.sided",paired = TRUE,exact = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT_rRNA.wilcox"), face="bold",fill="Specimen")
    } else if(m>0){
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "greater",paired = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT_rRNA.wilcox.greater"), face="bold",fill="Specimen")
    } else {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "less",paired = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT_rRNA.wilcox.test.less"), face="bold",fill="Batch")
    }
  }
  
  #comparison between CRC and NC in patients and healthy donor
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/10.RNA_biotypes/20210820_version2/")
    specimen_biotype <- read.csv("Specimen_biotypes.csv",header = TRUE, row.names = 1)
    
    #specimen_biotype <- specimen_biotype[c("CRC-PKU-27-PBMC","CRC-PKU-28-PBMC","CRC-PKU-29-PBMC","CRC-PKU-30-PBMC","CRC-PKU-32-PBMC",
    #                                       "CRC-PKU-27-pico","CRC-PKU-28-pico","CRC-PKU-29-pico","CRC-PKU-30-pico","CRC-PKU-32-pico",
    #                                       "CRC-PKU-27-T","CRC-PKU-28-T","CRC-PKU-29-T","CRC-PKU-30-T","CRC-PKU-32-T"),]
    
    specimen_biotype <- specimen_biotype[c("CRC-PKU-27-PBMC","CRC-PKU-28-PBMC","CRC-PKU-29-PBMC","CRC-PKU-30-PBMC","CRC-PKU-31-PBMC","CRC-PKU-32-PBMC",
                                           "NC-PKU-9-PBMC","NC-PKU-10-PBMC","NC-PKU-11-PBMC","NC-PKU-12-PBMC","NC-PKU-13-PBMC","NC-PKU-14-PBMC",
                                           "CRC-2280766",
                                           "CRC-2327118",
                                           "CRC-2334660",
                                           "CRC-2362297",
                                           "CRC-2372270",
                                           "CRC-2380493",
                                           "CRC-2382513",
                                           "CRC-2383018",
                                           "CRC-2383387",
                                           "CRC-2383420",
                                           "CRC-2383555",
                                           "CRC-2383730",
                                           "CRC-2384058",
                                           "CRC-2384076",
                                           "CRC-2384108",
                                           "CRC-2384226",
                                           "CRC-2384272",
                                           "CRC-2384377",
                                           "CRC-2384422",
                                           "CRC-2384541",
                                           "CRC-2384548",
                                           "CRC-2384796",
                                           "CRC-2386412",
                                           "CRC-2386612",
                                           "CRC-2387327",
                                           "CRC-2388859",
                                           "CRC-2388932",
                                           "CRC-2389212",
                                           "CRC-2390502",
                                           "CRC-2390539",
                                           "CRC-2390716",
                                           "CRC-2391097",
                                           "CRC-2393807",
                                           "CRC-2396217",
                                           "CRC-2396837",
                                           "CRC-2397180",
                                           "CRC-2398897",
                                           "CRC-2399129",
                                           "CRC-2405108",
                                           "CRC-2406345",
                                           "CRC-2406552",
                                           "CRC-2407821",
                                           "CRC-2408137",
                                           "CRC-2408162",
                                           "CRC-2408206",
                                           "CRC-2408346",
                                           "CRC-2408784",
                                           "CRC-2410759",
                                           "CRC-2410966",
                                           "CRC-2411883",
                                           "CRC-2411896",
                                           "CRC-2412048",
                                           "CRC-2412066",
                                           "CRC-9299429",
                                           "CRC-PKU-10-pico",
                                           "CRC-PKU-12-pico",
                                           "CRC-PKU-13-pico",
                                           "CRC-PKU-14-pico",
                                           "CRC-PKU-15-pico",
                                           "CRC-PKU-16-pico",
                                           "CRC-PKU-17-pico",
                                           "CRC-PKU-18-pico",
                                           "CRC-PKU-19-pico",
                                           "CRC-PKU-20-pico",
                                           "CRC-PKU-21-pico",
                                           "CRC-PKU-22-pico",
                                           "CRC-PKU-23-pico",
                                           "CRC-PKU-24-pico",
                                           "CRC-PKU-25-pico",
                                           "CRC-PKU-26-pico",
                                           "CRC-PKU-27-pico",
                                           "CRC-PKU-28-pico",
                                           "CRC-PKU-29-pico",
                                           "CRC-PKU-2-pico",
                                           "CRC-PKU-30-pico",
                                           "CRC-PKU-31-pico",
                                           "CRC-PKU-32-pico",
                                           "CRC-PKU-33-pico",
                                           "CRC-PKU-34-pico",
                                           "CRC-PKU-35-pico",
                                           "CRC-PKU-36-pico",
                                           "CRC-PKU-37-pico",
                                           "CRC-PKU-38-pico",
                                           "CRC-PKU-39-pico",
                                           "CRC-PKU-3-pico",
                                           "CRC-PKU-40-pico",
                                           "CRC-PKU-41-pico",
                                           "CRC-PKU-4-pico",
                                           "CRC-PKU-5-pico",
                                           "CRC-PKU-6-pico",
                                           "CRC-PKU-8-pico",
                                           "CRC-PKU-9-pico",
                                           "CRC-PKU-mix1-1-pico",
                                           "CRC-PKU-mix1-2-pico",
                                           "CRC-PKU-mix1-3-pico",
                                           "NC-PKU-10-pico",
                                           "NC-PKU-11-pico",
                                           "NC-PKU-12-pico",
                                           "NC-PKU-13-pico",
                                           "NC-PKU-14-pico",
                                           "NC-PKU-2359973",
                                           "NC-PKU-2392860",
                                           "NC-PKU-2392931",
                                           "NC-PKU-2392932",
                                           #"NC-PKU-2392993",
                                           "NC-PKU-2394780",
                                           "NC-PKU-2396996",
                                           "NC-PKU-2397512",
                                           "NC-PKU-2397878",
                                           "NC-PKU-2399728",
                                           "NC-PKU-28",
                                           "NC-PKU-33",
                                           "NC-PKU-34",
                                           "NC-PKU-37",
                                           "NC-PKU-48",
                                           "NC-PKU-50",
                                           "NC-PKU-56",
                                           "NC-PKU-69",
                                           "NC-PKU-8111237",
                                           "NC-PKU-9-pico",
                                           "NC-PKU-mix1-1-pico",
                                           "NC-PKU-mix1-2-pico",
                                           "NC-PKU-mix1-3-pico",
                                           "NC-PKU-mix15-pico",
                                           "NC-PKU-mix16-pico",
                                           "NC-PKU-mix17-pico",
                                           "NC-PKU-mix18-pico",
                                           "NC-PKU-mix19-pico",
                                           
                                           "NC-PKU-mix20-pico",
                                           "NC-PKU-mix2-1-pico",
                                           "NC-PKU-mix21-pico",
                                           "NC-PKU-mix2-2-pico",
                                           "NC-PKU-mix22-pico",
                                           "NC-PKU-mix24-pico",
                                           "NC-PKU-mix25-pico",
                                           "NC-PKU-mix26-pico",
                                           "NC-PKU-mix27-pico",
                                           "NC-PKU-mix28-pico",
                                           "NC-PKU-mix29-pico",
                                           
                                           "NC-PKU-mix30-pico",
                                           "NC-PKU-mix3-1-pico",
                                           "NC-PKU-mix3-2-pico",
                                           "NC-PKU-mix3-3-pico",
                                           "NC-PKU-mix3-4-pico",
                                           "NC-PKU-mix3-5-pico",
                                           
                                           "NC-PKU-mix4-1-pico",
                                           "NC-PKU-mix4-2-pico",
                                           
                                           "NC-PKU-mix5-1-pico",
                                           "NC-PKU-mix5-2-pico",
                                           "NC-PKU-mix5-3-pico",
                                           
                                           "NC-PKU-mix6-1-pico",
                                           "NC-PKU-mix6-2-pico",
                                           "NC-PKU-mix6-3-pico",
                                           "NC-PKU-mix6-4-pico",
                                           "NC-PKU-mix6-5-pico",
                                           "NC-PKU-mix6-6-pico",
                                           
                                           "NC-PKU-mix7-1-pico",
                                           "NC-PKU-mix7-2-pico",
                                           
                                           "NC-PKU-mix8-1-pico",
                                           "NC-PKU-mix8-2-pico",
                                           "NC-PKU-mix8-3-pico",
                                           "NC-PKU-mix8-4-pico",
                                           
                                           "NC-PKU-mix9-1-pico",
                                           "NC-PKU-mix9-2-pico",
                                           "NC-PKU-mix9-pico",
                                           
                                           #"CRC-PKU-27-pico","CRC-PKU-28-pico","CRC-PKU-29-pico","CRC-PKU-30-pico","CRC-PKU-31-pico","CRC-PKU-32-pico","CRC-PKU-34-pico","CRC-PKU-35-pico","CRC-PKU-36-pico","CRC-PKU-37-pico","CRC-PKU-38-pico","CRC-PKU-39-pico","CRC-PKU-40-pico","CRC-PKU-41-pico",
                                           #"NC-PKU-9-pico","NC-PKU-10-pico","NC-PKU-11-pico","NC-PKU-12-pico","NC-PKU-13-pico","NC-PKU-14-pico",
                                           "CRC-PKU-27-T","CRC-PKU-28-T","CRC-PKU-29-T","CRC-PKU-30-T","CRC-PKU-32-T","CRC-PKU-34-T","CRC-PKU-35-T","CRC-PKU-36-T","CRC-PKU-37-T","CRC-PKU-38-T","CRC-PKU-39-T","CRC-PKU-40-T","CRC-PKU-41-T",
                                           "CRC-PKU-27-N","CRC-PKU-28-N","CRC-PKU-29-N","CRC-PKU-30-N","CRC-PKU-32-N","CRC-PKU-34-N","CRC-PKU-35-N","CRC-PKU-36-N","CRC-PKU-37-N","CRC-PKU-38-N","CRC-PKU-39-N","CRC-PKU-40-N","CRC-PKU-41-N"),]
    
    specimen_biotype$MT <- specimen_biotype$MT_tRNA+specimen_biotype$MT_mRNA+specimen_biotype$MT_rRNA+specimen_biotype$MT_lncRNA+specimen_biotype$MT_exon
    specimen_biotype$sncRNA <- specimen_biotype$tRNA+specimen_biotype$srpRNA+specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA
    specimen_biotype$cleaned <- specimen_biotype$genome+specimen_biotype$circRNA+specimen_biotype$unmapped
    
    specimen_biotype$HG <- specimen_biotype$circRNA+specimen_biotype$genome-specimen_biotype$MT
    specimen_biotype$unassigned <- specimen_biotype$unmapped-specimen_biotype$microbe
    
    specimen_biotype_all <- specimen_biotype[,c("MT","HG","microbe","unassigned")]
    specimen_biotype_MT <- specimen_biotype[,c("MT_tRNA","MT_mRNA","MT_rRNA","MT_lncRNA")]
    specimen_biotype_HG <- specimen_biotype[,c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","tucpRNA","misc")]
    specimen_biotype_sncRNA <- specimen_biotype[,c("tRNA","srpRNA","snoRNA","snRNA","Y_RNA")]
    
    specimen_biotype_all_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_all),margin = 1))
    specimen_biotype_MT_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_MT),margin = 1))
    specimen_biotype_HG_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_HG),margin = 1))
    specimen_biotype_sncRNA_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_sncRNA),margin = 1))
    
    #all
    specimen_biotype_forplot <- melt(specimen_biotype_all_ratio)
    specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_all_ratio),ncol(specimen_biotype_all_ratio))
    specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_all_ratio))
    specimen_biotype_forplot$ID <- as.character(lapply(strsplit(specimen_biotype_forplot$Sample,"-"),function(x) paste(x[1],x[2],x[3], sep ="-", collapse =NULL)))
    specimen_biotype_forplot$ID <- paste(specimen_biotype_forplot$ID,specimen_biotype_forplot$Type,sep ="-", collapse =NULL)
    specimen_biotype_forplot$ID <- as.character(lapply(strsplit(specimen_biotype_forplot$ID,"-"),function(x) paste(x[1],x[2],x[3],x[4], sep ="-", collapse =NULL)))
    colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type","ID")
    
    #all
    specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                               levels = c("MT","HG","microbe","unassigned"))
    #HG
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("circRNA","lncRNA","pseudogene","sncRNA","exon","mRNA"))
    #sncRNA
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("Y_RNA","tRNA","srpRNA","snoRNA","snRNA","tucpRNA"))
    #MT
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("MT_tRNA","MT_rRNA","MT_lncRNA","MT_exon","MT_mRNA"))
    
    specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-CRC","PBMC-Healthy",
                                                                                   "Plasma-CRC","Plasma-Healthy","Plasma-STAD",
                                                                                   "Tissue-CRC","Tissue-Healthy"))
    my_comparisons <- list(c("PBMC-CRC","PBMC-Healthy"),c("Plasma-CRC","Plasma-Healthy"),c("Tissue-CRC","Tissue-Healthy"))
    p <- ggplot(specimen_biotype_forplot[which(specimen_biotype_forplot$Biotype=="MT"),],aes(x=Type,y=Ratio,fill=Type))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("red","blue","red","blue","red","blue")) +
      geom_vline(xintercept = c(2.5,4.5),linetype = "dashed",color="grey")+
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
    
    m=1
    if(m==0) {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "two.sided",paired = TRUE,exact = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT.wilcox"), face="bold",fill="Specimen")
    } else if(m>0) {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "greater",paired = TRUE,exact = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT.wilcox.greater"), face="bold",fill="Specimen")
    } else {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "less",paired = TRUE,exact = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT.wilcox.test.less"), face="bold",fill="Specimen")
    }
  }
  
}

#RNA_biotype based on gencode v3
{
  #CRC
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/10.RNA_biotypes/20210824_version3/")
    specimen_biotype <- read.csv("Specimen_biotypes.csv",header = TRUE, row.names = 1)
    
    specimen_biotype <- specimen_biotype[c("CRC-PKU-27-PBMC","CRC-PKU-28-PBMC","CRC-PKU-29-PBMC","CRC-PKU-30-PBMC","CRC-PKU-32-PBMC",
                                           "CRC-PKU-27-pico","CRC-PKU-28-pico","CRC-PKU-29-pico","CRC-PKU-30-pico","CRC-PKU-32-pico",
                                           "CRC-PKU-27-T","CRC-PKU-28-T","CRC-PKU-29-T","CRC-PKU-30-T","CRC-PKU-32-T"),]
    
    specimen_biotype$MT <- specimen_biotype$MT_tRNA+specimen_biotype$MT_mRNA+specimen_biotype$MT_lncRNA+specimen_biotype$MT_exon
    specimen_biotype$sncRNA <- specimen_biotype$tRNA+specimen_biotype$srpRNA+specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA
    specimen_biotype$cleaned <- specimen_biotype$genome+specimen_biotype$circRNA+specimen_biotype$unmapped-specimen_biotype$MT_rRNA
    
    specimen_biotype$HG <- specimen_biotype$circRNA+specimen_biotype$genome-specimen_biotype$MT-specimen_biotype$MT_rRNA
    specimen_biotype$unassigned <- specimen_biotype$unmapped-specimen_biotype$microbe
    
    specimen_biotype_all <- specimen_biotype[,c("MT","HG","microbe","unassigned")]
    specimen_biotype_MT <- specimen_biotype[,c("MT_tRNA","MT_mRNA","MT_lncRNA")]
    specimen_biotype_HG <- specimen_biotype[,c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","tucpRNA","misc")]
    specimen_biotype_sncRNA <- specimen_biotype[,c("tRNA","srpRNA","snoRNA","snRNA","Y_RNA")]
    
    specimen_biotype_all_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_all),margin = 1))
    specimen_biotype_MT_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_MT),margin = 1))
    specimen_biotype_HG_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_HG),margin = 1))
    specimen_biotype_sncRNA_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_sncRNA),margin = 1))
    write.csv(specimen_biotype_all_ratio,"all_ratio.csv")
    
    #MT+HG+microbe+unassiganed
    {
      specimen_biotype_forplot <- melt(specimen_biotype_all)
      specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_all),ncol(specimen_biotype_all))
      specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_all))
      colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
      
      specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                                 levels = c("MT","HG","microbe","unassigned"))
      specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                     "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                     "Tissue-Healthy","Tissue-CRC"))
      #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
      all <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
        geom_bar(stat = "identity",position = "fill")+
        xlab("")+
        ylab("")+
        theme_bw()+
        theme(
          #legend.position="none",
          plot.margin = unit(x=c(10,10,10,10),units="pt"),
          legend.position="top",
          legend.direction="vertical",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.title = element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          #axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))+
        scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
        scale_fill_manual(values = c("#BE2625","#A0522D","#EE9A00","#006633","#8C8C8C","#5D478B"))
    }
    #MT
    {
      specimen_biotype_forplot <- melt(specimen_biotype_MT)
      specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_MT),ncol(specimen_biotype_MT))
      specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_MT))
      colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
      
      specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                                 levels = c("MT_tRNA","MT_rRNA","MT_lncRNA","MT_exon","MT_mRNA"))
      specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                     "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                     "Tissue-Healthy","Tissue-CRC"))
      #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
      MT <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
        geom_bar(stat = "identity",position = "fill")+
        xlab("")+
        ylab("")+
        theme_bw()+
        theme(
          #legend.position="none",
          plot.margin = unit(x=c(10,10,10,10),units="pt"),
          legend.position="top",
          legend.direction="vertical",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.title = element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          #axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))+
        scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
        scale_fill_npg()
        #scale_fill_manual(values = c("#B22222","#CC3232","#EE3B3B","#FF6666","#FFCCCC"))
    }
    #HG
    {
      specimen_biotype_forplot <- melt(specimen_biotype_HG)
      specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_HG),ncol(specimen_biotype_HG))
      specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_HG))
      colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
      
      specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                                 levels = c("circRNA","lncRNA","pseudogene","sncRNA","tucpRNA","misc","mRNA"))
      specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                     "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                     "Tissue-Healthy","Tissue-CRC"))
      specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
      HG <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
        geom_bar(stat = "identity",position = "fill")+
        xlab("")+
        ylab("")+
        theme_bw()+
        theme(
          #legend.position="none",
          plot.margin = unit(x=c(10,10,10,10),units="pt"),
          legend.position="top",
          legend.direction="vertical",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.title = element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          #axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))+
        scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
        #scale_fill_manual(values = c("#5D478B","#8968CD","#AB82FF","#CCCCFF","#D9D9F3","#F8F8FF"))
        #scale_fill_aaas()
        scale_fill_manual(values=c("#631779","#ee2200","#008b45","#3b4992","#008280","#bb1921","#5f559b"))
      
    }
    #sncRNA
    {
      specimen_biotype_forplot <- melt(specimen_biotype_sncRNA)
      specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_sncRNA),ncol(specimen_biotype_sncRNA))
      specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_sncRNA))
      colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type")
      
      specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                                 levels = c("Y_RNA","tRNA","srpRNA","snoRNA","snRNA"))
      specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                     "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                     "Tissue-Healthy","Tissue-CRC"))
      #specimen_biotype_forplot$Sample <- factor(specimen_biotype_forplot$Sample,levels=unique(specimen_biotype_forplot$Sample))
      sncRNA <- ggplot(specimen_biotype_forplot,aes(x=Type,y=Ratio,fill=Biotype))+
        geom_bar(stat = "identity",position = "fill")+
        xlab("")+
        ylab("")+
        theme_bw()+
        theme(
          #legend.position="none",
          plot.margin = unit(x=c(10,10,10,10),units="pt"),
          legend.position="top",
          legend.direction="vertical",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.title = element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          #axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))+
        scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01))+
        scale_fill_npg()
        #scale_fill_manual(values = c("#00688B","#65909A","#50A6C2","#009ACD","#ADD8E6","#F0F8FF"))
    }
    ggarrange(MT,all,HG,sncRNA,align = "h", ncol = 4, nrow =1)
  }
  #comparison between specimen in paired CRC
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/10.RNA_biotypes/20210820_version2/")
    specimen_biotype <- read.csv("Specimen_biotypes.csv",header = TRUE, row.names = 1)
    
    specimen_biotype <- specimen_biotype[c("CRC-PKU-27-PBMC","CRC-PKU-28-PBMC","CRC-PKU-29-PBMC","CRC-PKU-30-PBMC","CRC-PKU-32-PBMC",
                                           "CRC-PKU-27-pico","CRC-PKU-28-pico","CRC-PKU-29-pico","CRC-PKU-30-pico","CRC-PKU-32-pico",
                                           "CRC-PKU-27-T","CRC-PKU-28-T","CRC-PKU-29-T","CRC-PKU-30-T","CRC-PKU-32-T"),]
    
    specimen_biotype$MT <- specimen_biotype$MT_tRNA+specimen_biotype$MT_mRNA+specimen_biotype$MT_rRNA+specimen_biotype$MT_lncRNA+specimen_biotype$MT_exon
    specimen_biotype$sncRNA <- specimen_biotype$tRNA+specimen_biotype$srpRNA+specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA
    specimen_biotype$cleaned <- specimen_biotype$genome+specimen_biotype$circRNA+specimen_biotype$unmapped
    
    specimen_biotype$HG <- specimen_biotype$circRNA+specimen_biotype$genome-specimen_biotype$MT
    specimen_biotype$unassigned <- specimen_biotype$unmapped-specimen_biotype$microbe
    
    specimen_biotype_all <- specimen_biotype[,c("MT","HG","microbe","unassigned")]
    specimen_biotype_MT <- specimen_biotype[,c("MT_tRNA","MT_mRNA","MT_rRNA","MT_lncRNA")]
    specimen_biotype_HG <- specimen_biotype[,c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","tucpRNA","misc")]
    specimen_biotype_sncRNA <- specimen_biotype[,c("tRNA","srpRNA","snoRNA","snRNA","Y_RNA")]
    
    specimen_biotype_all_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_all),margin = 1))
    specimen_biotype_MT_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_MT),margin = 1))
    specimen_biotype_HG_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_HG),margin = 1))
    specimen_biotype_sncRNA_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_sncRNA),margin = 1))
    
    #MT
    specimen_biotype_forplot <- melt(specimen_biotype_MT_ratio)
    specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_MT_ratio),ncol(specimen_biotype_MT_ratio))
    specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_MT_ratio))
    specimen_biotype_forplot$ID <- as.character(lapply(strsplit(specimen_biotype_forplot$Sample,"-"),function(x) paste(x[1],x[2],x[3], sep ="-", collapse =NULL)))
    colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type","ID")
    
    #all
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("MT","HG","microbe","unassigned"))
    #HG
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("circRNA","lncRNA","pseudogene","sncRNA","tucpRNA","mRNA"))
    #sncRNA
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("Y_RNA","tRNA","srpRNA","snoRNA","snRNA"))
    #MT
    specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                               levels = c("MT_tRNA","MT_rRNA","MT_lncRNA","MT_exon","MT_mRNA"))
    
    specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-Healthy","PBMC-CRC",
                                                                                   "Plasma-Healthy","Plasma-CRC","Plasma-STAD",
                                                                                   "Tissue-Healthy","Tissue-CRC"))
    my_comparisons <- list(c("Plasma-CRC","PBMC-CRC"),c("Plasma-CRC","Tissue-CRC"))
    p <- ggplot(specimen_biotype_forplot[which(specimen_biotype_forplot$Biotype=="MT_rRNA"),],aes(x=Type,y=Ratio,fill=Type))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
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
    
    m=0
    if(m==0) {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "two.sided",paired = TRUE,exact = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT_rRNA.wilcox"), face="bold",fill="Specimen")
    } else if(m>0){
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "greater",paired = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT_rRNA.wilcox.greater"), face="bold",fill="Specimen")
    } else {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "less",paired = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT_rRNA.wilcox.test.less"), face="bold",fill="Batch")
    }
  }
  
  #comparison between CRC and NC in patients and healthy donor
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/10.RNA_biotypes/20210820_version2/")
    specimen_biotype <- read.csv("Specimen_biotypes.csv",header = TRUE, row.names = 1)
    
    #specimen_biotype <- specimen_biotype[c("CRC-PKU-27-PBMC","CRC-PKU-28-PBMC","CRC-PKU-29-PBMC","CRC-PKU-30-PBMC","CRC-PKU-32-PBMC",
    #                                       "CRC-PKU-27-pico","CRC-PKU-28-pico","CRC-PKU-29-pico","CRC-PKU-30-pico","CRC-PKU-32-pico",
    #                                       "CRC-PKU-27-T","CRC-PKU-28-T","CRC-PKU-29-T","CRC-PKU-30-T","CRC-PKU-32-T"),]
    
    specimen_biotype <- specimen_biotype[c("CRC-PKU-27-PBMC","CRC-PKU-28-PBMC","CRC-PKU-29-PBMC","CRC-PKU-30-PBMC","CRC-PKU-31-PBMC","CRC-PKU-32-PBMC",
                                           "NC-PKU-9-PBMC","NC-PKU-10-PBMC","NC-PKU-11-PBMC","NC-PKU-12-PBMC","NC-PKU-13-PBMC","NC-PKU-14-PBMC",
                                           "CRC-2280766",
                                           "CRC-2327118",
                                           "CRC-2334660",
                                           "CRC-2362297",
                                           "CRC-2372270",
                                           "CRC-2380493",
                                           "CRC-2382513",
                                           "CRC-2383018",
                                           "CRC-2383387",
                                           "CRC-2383420",
                                           "CRC-2383555",
                                           "CRC-2383730",
                                           "CRC-2384058",
                                           "CRC-2384076",
                                           "CRC-2384108",
                                           "CRC-2384226",
                                           "CRC-2384272",
                                           "CRC-2384377",
                                           "CRC-2384422",
                                           "CRC-2384541",
                                           "CRC-2384548",
                                           "CRC-2384796",
                                           "CRC-2386412",
                                           "CRC-2386612",
                                           "CRC-2387327",
                                           "CRC-2388859",
                                           "CRC-2388932",
                                           "CRC-2389212",
                                           "CRC-2390502",
                                           "CRC-2390539",
                                           "CRC-2390716",
                                           "CRC-2391097",
                                           "CRC-2393807",
                                           "CRC-2396217",
                                           "CRC-2396837",
                                           "CRC-2397180",
                                           "CRC-2398897",
                                           "CRC-2399129",
                                           "CRC-2405108",
                                           "CRC-2406345",
                                           "CRC-2406552",
                                           "CRC-2407821",
                                           "CRC-2408137",
                                           "CRC-2408162",
                                           "CRC-2408206",
                                           "CRC-2408346",
                                           "CRC-2408784",
                                           "CRC-2410759",
                                           "CRC-2410966",
                                           "CRC-2411883",
                                           "CRC-2411896",
                                           "CRC-2412048",
                                           "CRC-2412066",
                                           "CRC-9299429",
                                           "CRC-PKU-10-pico",
                                           "CRC-PKU-12-pico",
                                           "CRC-PKU-13-pico",
                                           "CRC-PKU-14-pico",
                                           "CRC-PKU-15-pico",
                                           "CRC-PKU-16-pico",
                                           "CRC-PKU-17-pico",
                                           "CRC-PKU-18-pico",
                                           "CRC-PKU-19-pico",
                                           "CRC-PKU-20-pico",
                                           "CRC-PKU-21-pico",
                                           "CRC-PKU-22-pico",
                                           "CRC-PKU-23-pico",
                                           "CRC-PKU-24-pico",
                                           "CRC-PKU-25-pico",
                                           "CRC-PKU-26-pico",
                                           "CRC-PKU-27-pico",
                                           "CRC-PKU-28-pico",
                                           "CRC-PKU-29-pico",
                                           "CRC-PKU-2-pico",
                                           "CRC-PKU-30-pico",
                                           "CRC-PKU-31-pico",
                                           "CRC-PKU-32-pico",
                                           "CRC-PKU-33-pico",
                                           "CRC-PKU-34-pico",
                                           "CRC-PKU-35-pico",
                                           "CRC-PKU-36-pico",
                                           "CRC-PKU-37-pico",
                                           "CRC-PKU-38-pico",
                                           "CRC-PKU-39-pico",
                                           "CRC-PKU-3-pico",
                                           "CRC-PKU-40-pico",
                                           "CRC-PKU-41-pico",
                                           "CRC-PKU-4-pico",
                                           "CRC-PKU-5-pico",
                                           "CRC-PKU-6-pico",
                                           "CRC-PKU-8-pico",
                                           "CRC-PKU-9-pico",
                                           "CRC-PKU-mix1-1-pico",
                                           "CRC-PKU-mix1-2-pico",
                                           "CRC-PKU-mix1-3-pico",
                                           "NC-PKU-10-pico",
                                           "NC-PKU-11-pico",
                                           "NC-PKU-12-pico",
                                           "NC-PKU-13-pico",
                                           "NC-PKU-14-pico",
                                           "NC-PKU-2359973",
                                           "NC-PKU-2392860",
                                           "NC-PKU-2392931",
                                           "NC-PKU-2392932",
                                           #"NC-PKU-2392993",
                                           "NC-PKU-2394780",
                                           "NC-PKU-2396996",
                                           "NC-PKU-2397512",
                                           "NC-PKU-2397878",
                                           "NC-PKU-2399728",
                                           "NC-PKU-28",
                                           "NC-PKU-33",
                                           "NC-PKU-34",
                                           "NC-PKU-37",
                                           "NC-PKU-48",
                                           "NC-PKU-50",
                                           "NC-PKU-56",
                                           "NC-PKU-69",
                                           "NC-PKU-8111237",
                                           "NC-PKU-9-pico",
                                           "NC-PKU-mix1-1-pico",
                                           "NC-PKU-mix1-2-pico",
                                           "NC-PKU-mix1-3-pico",
                                           "NC-PKU-mix15-pico",
                                           "NC-PKU-mix16-pico",
                                           "NC-PKU-mix17-pico",
                                           "NC-PKU-mix18-pico",
                                           "NC-PKU-mix19-pico",
                                           
                                           "NC-PKU-mix20-pico",
                                           "NC-PKU-mix2-1-pico",
                                           "NC-PKU-mix21-pico",
                                           "NC-PKU-mix2-2-pico",
                                           "NC-PKU-mix22-pico",
                                           "NC-PKU-mix24-pico",
                                           "NC-PKU-mix25-pico",
                                           "NC-PKU-mix26-pico",
                                           "NC-PKU-mix27-pico",
                                           "NC-PKU-mix28-pico",
                                           "NC-PKU-mix29-pico",
                                           
                                           "NC-PKU-mix30-pico",
                                           "NC-PKU-mix3-1-pico",
                                           "NC-PKU-mix3-2-pico",
                                           "NC-PKU-mix3-3-pico",
                                           "NC-PKU-mix3-4-pico",
                                           "NC-PKU-mix3-5-pico",
                                           
                                           "NC-PKU-mix4-1-pico",
                                           "NC-PKU-mix4-2-pico",
                                           
                                           "NC-PKU-mix5-1-pico",
                                           "NC-PKU-mix5-2-pico",
                                           "NC-PKU-mix5-3-pico",
                                           
                                           "NC-PKU-mix6-1-pico",
                                           "NC-PKU-mix6-2-pico",
                                           "NC-PKU-mix6-3-pico",
                                           "NC-PKU-mix6-4-pico",
                                           "NC-PKU-mix6-5-pico",
                                           "NC-PKU-mix6-6-pico",
                                           
                                           "NC-PKU-mix7-1-pico",
                                           "NC-PKU-mix7-2-pico",
                                           
                                           "NC-PKU-mix8-1-pico",
                                           "NC-PKU-mix8-2-pico",
                                           "NC-PKU-mix8-3-pico",
                                           "NC-PKU-mix8-4-pico",
                                           
                                           "NC-PKU-mix9-1-pico",
                                           "NC-PKU-mix9-2-pico",
                                           "NC-PKU-mix9-pico",
                                           
                                           #"CRC-PKU-27-pico","CRC-PKU-28-pico","CRC-PKU-29-pico","CRC-PKU-30-pico","CRC-PKU-31-pico","CRC-PKU-32-pico","CRC-PKU-34-pico","CRC-PKU-35-pico","CRC-PKU-36-pico","CRC-PKU-37-pico","CRC-PKU-38-pico","CRC-PKU-39-pico","CRC-PKU-40-pico","CRC-PKU-41-pico",
                                           #"NC-PKU-9-pico","NC-PKU-10-pico","NC-PKU-11-pico","NC-PKU-12-pico","NC-PKU-13-pico","NC-PKU-14-pico",
                                           "CRC-PKU-27-T","CRC-PKU-28-T","CRC-PKU-29-T","CRC-PKU-30-T","CRC-PKU-32-T","CRC-PKU-34-T","CRC-PKU-35-T","CRC-PKU-36-T","CRC-PKU-37-T","CRC-PKU-38-T","CRC-PKU-39-T","CRC-PKU-40-T","CRC-PKU-41-T",
                                           "CRC-PKU-27-N","CRC-PKU-28-N","CRC-PKU-29-N","CRC-PKU-30-N","CRC-PKU-32-N","CRC-PKU-34-N","CRC-PKU-35-N","CRC-PKU-36-N","CRC-PKU-37-N","CRC-PKU-38-N","CRC-PKU-39-N","CRC-PKU-40-N","CRC-PKU-41-N"),]
    
    specimen_biotype$MT <- specimen_biotype$MT_tRNA+specimen_biotype$MT_mRNA+specimen_biotype$MT_rRNA+specimen_biotype$MT_lncRNA+specimen_biotype$MT_exon
    specimen_biotype$sncRNA <- specimen_biotype$tRNA+specimen_biotype$srpRNA+specimen_biotype$snoRNA+specimen_biotype$snRNA+specimen_biotype$Y_RNA
    specimen_biotype$cleaned <- specimen_biotype$genome+specimen_biotype$circRNA+specimen_biotype$unmapped
    
    specimen_biotype$HG <- specimen_biotype$circRNA+specimen_biotype$genome-specimen_biotype$MT
    specimen_biotype$unassigned <- specimen_biotype$unmapped-specimen_biotype$microbe
    
    specimen_biotype_all <- specimen_biotype[,c("MT","HG","microbe","unassigned")]
    specimen_biotype_MT <- specimen_biotype[,c("MT_tRNA","MT_mRNA","MT_rRNA","MT_lncRNA")]
    specimen_biotype_HG <- specimen_biotype[,c("mRNA","lncRNA","circRNA","pseudogene","sncRNA","tucpRNA","misc")]
    specimen_biotype_sncRNA <- specimen_biotype[,c("tRNA","srpRNA","snoRNA","snRNA","Y_RNA")]
    
    specimen_biotype_all_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_all),margin = 1))
    specimen_biotype_MT_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_MT),margin = 1))
    specimen_biotype_HG_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_HG),margin = 1))
    specimen_biotype_sncRNA_ratio <- as.data.frame(prop.table(as.matrix(specimen_biotype_sncRNA),margin = 1))
    
    #all
    specimen_biotype_forplot <- melt(specimen_biotype_all_ratio)
    specimen_biotype_forplot$Sample <- rep(rownames(specimen_biotype_all_ratio),ncol(specimen_biotype_all_ratio))
    specimen_biotype_forplot$Type <- rep(specimen_biotype$label,ncol(specimen_biotype_all_ratio))
    specimen_biotype_forplot$ID <- as.character(lapply(strsplit(specimen_biotype_forplot$Sample,"-"),function(x) paste(x[1],x[2],x[3], sep ="-", collapse =NULL)))
    specimen_biotype_forplot$ID <- paste(specimen_biotype_forplot$ID,specimen_biotype_forplot$Type,sep ="-", collapse =NULL)
    specimen_biotype_forplot$ID <- as.character(lapply(strsplit(specimen_biotype_forplot$ID,"-"),function(x) paste(x[1],x[2],x[3],x[4], sep ="-", collapse =NULL)))
    colnames(specimen_biotype_forplot) <- c("Biotype","Ratio","Sample","Type","ID")
    
    #all
    specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
                                               levels = c("MT","HG","microbe","unassigned"))
    #HG
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("circRNA","lncRNA","pseudogene","sncRNA","exon","mRNA"))
    #sncRNA
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("Y_RNA","tRNA","srpRNA","snoRNA","snRNA","tucpRNA"))
    #MT
    #specimen_biotype_forplot$Biotype <- factor(specimen_biotype_forplot$Biotype,
    #                                           levels = c("MT_tRNA","MT_rRNA","MT_lncRNA","MT_exon","MT_mRNA"))
    
    specimen_biotype_forplot$Type <- factor(specimen_biotype_forplot$Type,levels=c("PBMC-CRC","PBMC-Healthy",
                                                                                   "Plasma-CRC","Plasma-Healthy","Plasma-STAD",
                                                                                   "Tissue-CRC","Tissue-Healthy"))
    my_comparisons <- list(c("PBMC-CRC","PBMC-Healthy"),c("Plasma-CRC","Plasma-Healthy"),c("Tissue-CRC","Tissue-Healthy"))
    p <- ggplot(specimen_biotype_forplot[which(specimen_biotype_forplot$Biotype=="MT"),],aes(x=Type,y=Ratio,fill=Type))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("red","blue","red","blue","red","blue")) +
      geom_vline(xintercept = c(2.5,4.5),linetype = "dashed",color="grey")+
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
    
    m=1
    if(m==0) {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "two.sided",paired = TRUE,exact = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT.wilcox"), face="bold",fill="Specimen")
    } else if(m>0) {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "greater",paired = TRUE,exact = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT.wilcox.greater"), face="bold",fill="Specimen")
    } else {
      p <- p+geom_line(aes(group = ID))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "less",paired = TRUE,exact = TRUE),
                           label = "p.signif"
        )+labs(x="",y="Ratio (%)",title=paste0("MT.wilcox.test.less"), face="bold",fill="Specimen")
    }
  }
  
}

#Correlation between tissue and plasma at pathway level
{
  #single scatter plot
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/DE")
  DE_plasma <- read.csv("plasma/CRCvsNC_edger_glmlrt-plasma.txt",sep = "\t",header = T)
  DE_plasma$index <- rownames(DE_plasma)
  colnames(DE_plasma) <- c("log2FoldChange.plasma","logCPM","LR","pvalue","padj","baseMean","index")
  DE_tissue <- read.csv("tissue/TumorvsNormal_edger_glmlrt-tissue.txt",sep = "\t",header = T)
  DE_tissue$index <- rownames(DE_tissue)
  colnames(DE_tissue) <- c("log2FoldChange.tissue","logCPM","LR","pvalue","padj","baseMean","index")
  DE_PBMC <- read.csv("PBMC/CRCvsNC_edger_glmlrt-PBMC.txt",sep = "\t",header = T)
  DE_PBMC$index <- rownames(DE_PBMC)
  colnames(DE_PBMC) <- c("log2FoldChange.PBMC","logCPM","LR","pvalue","padj","baseMean","index")
  
  DE_plasma_tissue <- inner_join(DE_plasma,DE_tissue, by = c("index"="index") )
  DE_plasma_tissue_PBMC <- inner_join(DE_plasma_tissue,DE_PBMC, by = c("index"="index") )
  
  DE_plasma_tissue_PBMC$ENSG <- as.character(lapply(strsplit(DE_plasma_tissue_PBMC$index,".",fixed = TRUE),function(x) x[1]))
  
  library("ggpubr")
  ggscatter(DE_plasma_tissue_PBMC[RP_RNA$Gene.ID,], x = "log2FoldChange.tissue", y = "log2FoldChange.PBMC", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson")
  
  # cor and pvalue bubble plot
  library(corrplot)
  library(RColorBrewer)
  pathway_forcorr <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation/Correlation between plasma expression and tissue immune fraction.csv",header = TRUE, row.names = 1)
  M <-cor(pathway_forcorr)
  library("Hmisc")
  col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                 "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                 "#4393C3", "#2166AC", "#053061")))
  # Insignificant correlation are crossed
  res2 <- rcorr(as.matrix(pathway_forcorr))
  corrplot(corr = res2$r,col = col2(200),tl.col="black",type="lower", order="original",tl.pos = "ld",#tl.cex=0.7,tl.srt = 45,
           p.mat = res2$P, sig.level = 0.05,insig = "blank")  # insig = "blank"

  #multiple scatter plot and cor-pvalue matrix
  {
  library(pspearman)
  spearman_CI <- function(x, y, alpha = 0.05){
    rs <- cor(x, y, method = "spearman", use = "complete.obs")
    n <- sum(complete.cases(x, y))
    sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
  }
  
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation/Plasma and PBMC immune fraction_spearman/"
  cor_methods <- "spearman"
  forcor <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation/Correlation between plasma expression and PBMC immune fraction.csv",header = TRUE, row.names = 1)
  forcor <- forcor[-c(1,2),] #for Correlation between plasma expression and tissue immune fraction.csv
  
  #forcor <- forcor[,1:10]
  
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = ncol(forcor), clear = FALSE, width= 60)
  
  result_final <- as.data.frame(matrix(numeric(0),ncol=7))
  colnames(result_final) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
  
  result <- as.data.frame(matrix(numeric(0),ncol=7))
  colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
  i=1
  while(i<=(ncol(forcor)-1)){
    #print(i)
    j=1
    while(j<=(ncol(forcor)-i)){
      
      if(cor_methods=="pearson"){
        r_twosided <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods)
        if (is.na(r_twosided$estimate)) {
          r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
        } else if (r_twosided$estimate>0){
          r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "greater")
        } else if (r_twosided$estimate<0) {
          r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "less")
        } else {
          r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
        }
        result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,r$conf.int[1],r$conf.int[2],r_twosided$p.value)
      } else if(cor_methods=="spearman") {
        r_twosided <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
        if (is.na(r_twosided$estimate)) {
          r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
        } else if(r_twosided$estimate>0){
          r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "greater", approximation = "exact")
        } else if (r_twosided$estimate<0) {
          r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "less", approximation = "exact")
        } else {
          r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
        }
        result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,spearman_CI(forcor[,i],forcor[,i+j])[1],spearman_CI(forcor[,i],forcor[,i+j])[2],r_twosided$p.value)
      }
      
      colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
      result <- rbind(result,result_tmp)
      j=j+1
    }
    if(i%%20==0){
      result_final <- rbind(result_final,result)
      write.csv(result,paste0(output_dir,i/20,"_Result20.csv"))
      result <- as.data.frame(matrix(numeric(0),ncol=7))
      colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
      i=i+1
    } else if(i==(ncol(forcor)-1)){
      result_final <- rbind(result_final,result)
      write.csv(result,paste0(output_dir,i%/%20+1,"_Result",i%%20,".csv"))
      i=i+1
    } else {
      i=i+1
    }
    pb$tick()
    Sys.sleep(1 / 100)
  }
  write.csv(result_final,paste0(output_dir,"Result_final.csv"))
  }
  
  #circos plot for overall correlation
  install.packages("circlize")
  library(circlize)
  
  mat = matrix(1:9, 3)
  rownames(mat) = letters[1:3]
  colnames(mat) = LETTERS[1:3]
  mat
  
  df <- melt(mat)
  colnames(df) <- c("from","to","value")
  
  chordDiagram(mat)
  chordDiagram(df,link.visible = df[[3]] >= 7)
  
  cor_df <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation/correlation_Spearman_forcircos.csv",header=TRUE)
  pathway_group <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation/pathway_group.csv",header = TRUE)
  #nm = unique(unlist(pathway_group$Pathway))
  #group = structure(pathway_group$group, names = as.character(nm))
  #track_color = structure(pathway_group$color, names = as.character(nm))
  
  nm = c("Plasma","PBMC","Normal","Tumor")
  group = structure(c("Plasma","PBMC","Normal","Tumor"), names = nm)
  color = structure(c("#8B6969","#FFFFFF","#3D59AB","#EE0000"), names = nm)
  
  cor_circoplot <- as.data.frame(matrix(numeric(0),ncol=4,nrow = nrow(cor_df)))
  colnames(cor_circoplot) <- c("from","to","value","group")
  
  cor_circoplot$from <- as.character(lapply(strsplit(as.character(cor_df$Comparison.group)," vs ",fixed = TRUE),function(x) x[1]))
  cor_circoplot$to <- as.character(lapply(strsplit(as.character(cor_df$Comparison.group)," vs ",fixed = TRUE),function(x) x[2]))
  cor_circoplot$value <- cor_df$significant.correlated
  #cor_circoplot$value <- gsub("-1","3",cor_circoplot$value)
  cor_circoplot$group <- cor_df$DataNames
  cor_circoplot[is.na(cor_circoplot)] <- 3
  chordDiagram(cor_circoplot,grid.col = color,annotationTrack = "grid",link.visible = cor_circoplot[[3]] == 1, link.lwd = 0.2, link.lty = 3, link.border = "#000000")
  
  #single scatter plot and cor-pvalue matrix
  forcor <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation/Correlation between plasma expression and PBMC immune fraction.csv",header = TRUE, row.names = 1)
  forcor <- forcor[-c(1,2),] #for Correlation between plasma expression and tissue immune fraction.csv
  ggscatter(forcor, x = "RP.mRNA.Plasma", y = "T.cells.CD4.memory.activated.ImmuneFraction", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "spearman", label.x = 7, label.y = 0.05,label.sep = "\n",size = 8))+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
}

#correlation between plasma and tissue at sample level
{
  counts_plasma <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/plasma/TPM_intron_spanning_noMTRNA.txt",sep="\t",header = T, row.names = 1)
  counts_tissue <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/tissue/tissue_TPM_intron_spanning_noMTRNA.txt",sep="\t",header = T, row.names = 1)
  counts_PBMC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/PBMC/TPM_intron_spanning_noMTRNA.txt",sep="\t",header = T, row.names = 1)
 
  { 
  plasma_forcor <- counts_plasma_tissue[,c("CRC.PKU.27.pico","CRC.PKU.28.pico","CRC.PKU.30.pico",
                                           "CRC.PKU.32.pico","CRC.PKU.35.pico","CRC.PKU.36.pico",
                                           "CRC.PKU.37.pico","CRC.PKU.38.pico","CRC.PKU.39.pico","CRC.PKU.40.pico")]
  colnames(plasma_forcor) <- gsub(".","-",fixed=TRUE,colnames(plasma_forcor))
  colnames(plasma_forcor) <- gsub("-pico","-plasma",fixed=TRUE,colnames(plasma_forcor))
  
  tumor_forcor <- counts_plasma_tissue[,c("CRC.PKU.27.T","CRC.PKU.28.T","CRC.PKU.30.T",
                                          "CRC.PKU.32.T","CRC.PKU.35.T","CRC.PKU.36.T",
                                          "CRC.PKU.37.T","CRC.PKU.38.T","CRC.PKU.39.T","CRC.PKU.40.T")]
  colnames(tumor_forcor) <- gsub(".","-",fixed=TRUE,colnames(tumor_forcor))
  colnames(tumor_forcor) <- gsub("-T","-tumor",fixed=TRUE,colnames(tumor_forcor))
 
  normal_forcor <- counts_plasma_tissue[,c("CRC.PKU.27.N","CRC.PKU.28.N","CRC.PKU.30.N",
                                           "CRC.PKU.32.N","CRC.PKU.35.N","CRC.PKU.36.N",
                                           "CRC.PKU.37.N","CRC.PKU.38.N","CRC.PKU.39.N","CRC.PKU.40.N")]
  colnames(normal_forcor) <- gsub(".","-",fixed=TRUE,colnames(normal_forcor))
  colnames(normal_forcor) <- gsub("-N","-normal",fixed=TRUE,colnames(normal_forcor))
  
  plasma_tissue_paired <- cbind(plasma_forcor,tumor_forcor,normal_forcor)
  plasma_tissue_paired_ensembl <- plasma_tissue_paired[grep("ENSG",rownames(plasma_tissue_paired)),]
  
  # cor and pvalue bubble plot
  library(corrplot)
  library(RColorBrewer)
  
  M <-cor(plasma_tissue_paired_ensembl)
  library("Hmisc")
  col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                 "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                 "#4393C3", "#2166AC", "#053061")))
  # Insignificant correlation are crossed
  res2 <- rcorr(as.matrix(plasma_tissue_paired_ensembl),type="pearson")
  corrplot(corr = res2$r,col = col2(200),tl.col="black",type="lower", order="original",tl.pos = "ld",#tl.cex=0.7,tl.srt = 45,
           p.mat = res2$P, sig.level = 0.0000001,insig = "blank")  # insig = "blank"
  }
  
  #plasma vs PBMC vs tumor
  {
  plasma_PBMC_CRC_forcor <- counts_plasma[,c("CRC.PKU.29.pico","CRC.PKU.27.pico","CRC.PKU.28.pico",
                                             "CRC.PKU.30.pico","CRC.PKU.32.pico")]
  colnames(plasma_PBMC_CRC_forcor) <- gsub(".","-",fixed=TRUE,colnames(plasma_PBMC_CRC_forcor))
  colnames(plasma_PBMC_CRC_forcor) <- gsub("-pico","-plasma",fixed=TRUE,colnames(plasma_PBMC_CRC_forcor))
  
  PBMC_CRC_forcor <- counts_PBMC[,c("CRC.PKU.29.PBMC","CRC.PKU.27.PBMC","CRC.PKU.28.PBMC",
                                    "CRC.PKU.30.PBMC","CRC.PKU.32.PBMC")]
  colnames(PBMC_CRC_forcor) <- gsub(".","-",fixed=TRUE,colnames(PBMC_CRC_forcor))
  
  tumor_forcor <- counts_plasma_tissue[,c("CRC.PKU.29.T","CRC.PKU.27.T","CRC.PKU.28.T","CRC.PKU.30.T",
                                          "CRC.PKU.32.T")]
  colnames(tumor_forcor) <- gsub(".","-",fixed=TRUE,colnames(tumor_forcor))
  colnames(tumor_forcor) <- gsub("-T","-tumor",fixed=TRUE,colnames(tumor_forcor))
  
  plasma_PBMC_tumor_paired <- cbind(plasma_PBMC_CRC_forcor,PBMC_CRC_forcor,tumor_forcor)
  plasma_PBMC_tumor_paired_ensembl <- plasma_PBMC_tumor_paired[grep("ENSG",rownames(plasma_PBMC_tumor_paired)),]
  
  # cor and pvalue bubble plot
  library(corrplot)
  library(RColorBrewer)
  
  M <-cor(plasma_PBMC_tumor_paired_ensembl)
  library("Hmisc")
  col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                 "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                 "#4393C3", "#2166AC", "#053061")))
  # Insignificant correlation are crossed
  res2 <- rcorr(as.matrix(plasma_PBMC_tumor_paired_ensembl),type="pearson")
  corrplot(corr = res2$r,col = col2(200),tl.col="black",tl.pos = "l",type="lower", order="original",diag = TRUE,#tl.cex=0.7,tl.srt = 45,tl.pos = "ld",type="lower"
           p.mat = res2$P, sig.level = 0.01,insig = "blank",cl.lim = c(0,1))  # insig = "blank"
  
  ggscatter(plasma_PBMC_tumor_paired_ensembl, x = "CRC-PKU-27-plasma", y = "CRC-PKU-28-plasma", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      plot.margin = unit(c(20,20,20,20),"pt"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
  
  ggscatter(plasma_PBMC_tumor_paired_ensembl, x = "CRC-PKU-27-PBMC", y = "CRC-PKU-27-plasma", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "pearson",label.sep = "\n",size = 8))+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      plot.margin = unit(c(20,20,20,20),"pt"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
  
  ggscatter(plasma_PBMC_tumor_paired_ensembl, x = "CRC-PKU-27-PBMC", y = "CRC-PKU-28-PBMC", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "pearson",label.sep = "\n",size = 8))+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      plot.margin = unit(c(20,20,20,20),"pt"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
  
  ggscatter(plasma_PBMC_tumor_paired_ensembl, x = "CRC-PKU-27-tumor", y = "CRC-PKU-27-plasma", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "pearson",label.sep = "\n",size = 8))+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      plot.margin = unit(c(20,20,20,20),"pt"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
  
  ggscatter(plasma_PBMC_tumor_paired_ensembl, x = "CRC-PKU-27-tumor", y = "CRC-PKU-28-tumor", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "pearson",label.sep = "\n",size = 8))+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      plot.margin = unit(c(20,20,20,20),"pt"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
  }
  
}

#differential expression between CRC and NC (Tumor and Normal) in plasma/PBMC/tissue, rank by foldchange genes (FDR<0.05) and then enrich at pathway level by GSEA
{
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/specimen_DE_GSEA/")
  mat_raw <- read.table("./matrix/PBMC_featurecount_intron_spanning_noMTRNA.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
  des <- read.csv("./group/PBMC.csv", header = TRUE, check.names=FALSE, sep=',')
  
  samples <- des$samples
  #samples <- gsub(".","-",samples,fixed=TRUE)
  #samples <- gsub("X","",samples,fixed=TRUE)
  group <- des$group
  batch <- des$batch
  
  i=1
  mat=as.data.frame(array(dim=c(length(rownames(mat_raw)),1)))
  while (i<=length(samples)) {
    temp <- mat_raw[,which(colnames(mat_raw)==samples[i])]
    #temp <- mat_raw[,grep(samples[i],colnames(mat_raw),fixed=TRUE)]
    #print(i)
    #print(which(colnames(mat_raw)==samples[i]))
    temp <- as.data.frame(temp)
    colnames(temp) <- samples[i]
    mat <- cbind(mat,temp)
    i=i+1
  }
  
  mat <- mat[,-1]
  rownames(mat) <- rownames(mat_raw)
  
  y <- DGEList(counts=mat, samples=samples, group=group)     #规定输入格式，必须包含这几个参数
  y <- calcNormFactors(y, method="TMM")                      #归一化处理
  y <- estimateDisp(y)
  test <- exactTest(y, pair = c("negative","positive"), dispersion = "auto")
  #https://support.bioconductor.org/p/64807/
  #广义线性模型可以针对多因素的情况，在这个场景下，只比较cancer与normal用exactTest就足够了
  {
    #design <- model.matrix(~0+group)   
    #y <- estimateDisp(y, design)          #广义线性模型计算离散度（common&trended&tagwise）
    #fit <- glmFit(y, design)              #差异表达分析函数                     
    #test <- glmLRT(fit, coef=2)           #差异表达分析函数 
  }
  
  
  
  res <- topTags(test, n=nrow(mat), sort.by='none')     #输出计算差异基因的结果
  res <- cbind(res$table, baseMean=2^(res$table$logCPM)) #输出差异基因结果加上baseMean一栏
  mapped_names <- colnames(res)                   #后面都是整理+输出
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  write.table(res, "./output/PBMC_edger_exact.csv", sep=',', quote=FALSE, row.names=TRUE) #输出文件，规定名字
}

#make ranked genelist
{
  res.sort <- res[sort(res$log2FoldChange,decreasing = T, index.return = T)$ix,]
  res.sort <- res.sort[grep("ENSG",rownames(res.sort)),]
  
  genelist <- as.numeric(res.sort$log2FoldChange)
  names(genelist) <- as.character(lapply(strsplit(rownames(res.sort),"\\|"), function(x) x[1]))
  names(genelist) <- as.character(lapply(strsplit(names(genelist),".",fixed = TRUE), function(x) x[1]))
  
}

#make geneset from KEGG
{
  PATH_ID_NAME <- read.csv("pathway/PATH_ID_NAME_KEGGplusHallmark.txt",header = TRUE, row.names = 1,sep = "\t",check.names=FALSE)
  #gmt <- PATH_ID_NAME[,c(4,1)] # local_modified
  gmt <- PATH_ID_NAME[,c(6,1)] # local_plusHallmarker
  colnames(gmt) <- c("ont","gene")
}

#GSEA by clusterprofile and plot
{
  egmt2 <- GSEA(genelist,minGSSize = 0, maxGSSize = 1000,
                pvalueCutoff = 1, TERM2GENE = gmt)
  write.table(egmt2@result,"./output/PBMC_GSEA.txt",quote = FALSE,sep = "\t")
  
  #gsea <- gseaplot2(egmt2,1)
  #pdf("Ribosome_PBMC_CRCvsNC.pdf")
  #plot(gsea)
  #dev.off()
}

#bubble plot
{
plasma <- read.csv("./output/plasma_all_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
plasma$Specimen <- "Plasma"
tissue <- read.csv("./output/tissue_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
tissue$Specimen <- "Tissue"
PBMC <- read.csv("./output/PBMC_GSEA.txt",header = TRUE,row.names = 1,sep = "\t")
PBMC$Specimen <- "PBMC"

forbubbleplot <- rbind(plasma,tissue,PBMC)

GSEA_res <- PBMC
filtered_up <- GSEA_res[(GSEA_res$p.adjust < 0.1) & (GSEA_res$NES > 0), ]
up <- head(filtered_up[order(filtered_up$NES,decreasing = TRUE),]$Description,100)
up <- factor(up,levels=up)
filtered_down <- GSEA_res[(GSEA_res$p.adjust < 0.1) & (GSEA_res$NES < 0), ]
down <- head(filtered_down[order(filtered_down$NES,decreasing = FALSE),]$Description,10)
down <- factor(down,levels=down)
top20 <- fct_c(up,down)

genes <- top20
genes <- c(
  "Metabolics","Glycolysis Gluconeogenesis","Lipoprotein metabolism","Biosynthesis of amino acids","Protein export","Protein digestion absorption",
  "Pro-inflamatory cytokines&chemokines","Inflammatory response","Unfolded protein response","Primary immunodeficiency","Human immunodeficiency virus 1 infection",
  "Epithelial mesenchymal transition","G2M checkpoint","Reactive oxigen species","Hypoxia","Angiogensis",
  "Apoptosis","Resisting Cell Death","Sustaining Proliferative Signaling","Extension of telomeres","Protein secretion",
  "MT-RP-RNA","RP-mRNA","Ribosome"
)


forbubbleplot_plot <- filter(forbubbleplot, Description %in% genes)
forbubbleplot_plot[forbubbleplot_plot$p.adjust>=0.1,"p.adjust"] <- NA

forbubbleplot_plot$Description <- factor(forbubbleplot_plot$Description,levels = genes)
forbubbleplot_plot$Specimen <- factor(forbubbleplot_plot$Specimen,levels = c("Tissue","PBMC","Plasma"))
ggplot(forbubbleplot_plot,aes(x=forbubbleplot_plot$Specimen,y=forbubbleplot_plot$Description))+
  geom_point(aes(size=-1*log10(forbubbleplot_plot$p.adjust),color=forbubbleplot_plot$NES))+
  scale_colour_gradient(low="blue",high="red")+
  labs(color="NES",
       size=expression(-log[10](p.adjust)),
       x="Specimen")+
  theme_bw()+
  theme(axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
        axis.text.x = element_text(size=rel(1.3),face="bold",colour = "black"),
        axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
        axis.title.y = element_blank())
}
}

#differential expression between Tumor(cancer patients) and Normal(healthy donors) in public datasets, select genes (FDR<0.05) and then enrich at pathway level by KEGG
{
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/09.Published_data/")
    mat_raw <- read.table("./matrix/PBMC_featurecount_intron_spanning_noMTRNA.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    des <- read.csv("./group/PBMC.csv", header = TRUE, check.names=FALSE, sep=',')
    
    samples <- des$samples
    #samples <- gsub(".","-",samples,fixed=TRUE)
    #samples <- gsub("X","",samples,fixed=TRUE)
    group <- des$group
    batch <- des$batch
    
    i=1
    mat=as.data.frame(array(dim=c(length(rownames(mat_raw)),1)))
    while (i<=length(samples)) {
      temp <- mat_raw[,which(colnames(mat_raw)==samples[i])]
      #temp <- mat_raw[,grep(samples[i],colnames(mat_raw),fixed=TRUE)]
      #print(i)
      #print(which(colnames(mat_raw)==samples[i]))
      temp <- as.data.frame(temp)
      colnames(temp) <- samples[i]
      mat <- cbind(mat,temp)
      i=i+1
    }
    
    mat <- mat[,-1]
    rownames(mat) <- rownames(mat_raw)
    
    y <- DGEList(counts=mat, samples=samples, group=group)     #规定输入格式，必须包含这几个参数
    y <- calcNormFactors(y, method="TMM")                      #归一化处理
    y <- estimateDisp(y)
    test <- exactTest(y, pair = c("negative","positive"), dispersion = "auto")
    #https://support.bioconductor.org/p/64807/
    #广义线性模型可以针对多因素的情况，在这个场景下，只比较cancer与normal用exactTest就足够了
    {
      #design <- model.matrix(~0+group)   
      #y <- estimateDisp(y, design)          #广义线性模型计算离散度（common&trended&tagwise）
      #fit <- glmFit(y, design)              #差异表达分析函数                     
      #test <- glmLRT(fit, coef=2)           #差异表达分析函数 
    }
    
    
    
    res <- topTags(test, n=nrow(mat), sort.by='none')     #输出计算差异基因的结果
    res <- cbind(res$table, baseMean=2^(res$table$logCPM)) #输出差异基因结果加上baseMean一栏
    mapped_names <- colnames(res)                   #后面都是整理+输出
    for(i in 1:ncol(res)){
      if(colnames(res)[i] == 'logFC'){
        mapped_names[i] <- 'log2FoldChange'
      }else if(colnames(res)[i] == 'PValue'){
        mapped_names[i] <- 'pvalue'
      }else if(colnames(res)[i] == 'FDR') {
        mapped_names[i] <- 'padj'
      }else{
        mapped_names[i] <- colnames(res)[i]
      }
    }
    colnames(res) <- mapped_names
    write.table(res, "./output/PBMC_edger_exact.csv", sep=',', quote=FALSE, row.names=TRUE) #输出文件，规定名字
  }
  path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/09.Published_data/output/"
  docs <- dir(path)
  docs <- grep("_PrimaryTumor_vs_TissueNormal_edger_exact.csv",docs,value = TRUE)
  dataset_name <- as.character(lapply(strsplit(docs,"_edger_exact"),function(x) x[1]))
  print(dataset_name)
  i=1
  while(i<=length(dataset_name)){
    #select gene FDR < 0.05
    DE <- read.csv(paste0("output/",dataset_name[i],"_edger_exact.csv"),header = TRUE, row.names = 1)
    
    all_gene <- row.names(DE)
    all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
    filtered <- rownames(DE[DE$padj<0.05,])
    filtered <- as.character(lapply(strsplit(filtered,"\\."),function(x) x[1]))
    
    forenrich <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                       filters = "ensembl_gene_id",
                       values=filtered, mart= mart,useCache = FALSE)
    background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                        filters = "ensembl_gene_id",
                        values=all_gene, mart= mart,useCache = FALSE)
    
    forenrich <- na.omit(forenrich)
    background <- na.omit(background)
    
    KEGG_res <- enrichKEGG(
      forenrich$entrezgene_id,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      universe=as.character(background$entrezgene_id),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1,
      use_internal_data = FALSE
    )
    
    #dotplot(KEGG_res)
    
    write.csv(KEGG_res,paste0("output/",dataset_name[i],"_KEGG.csv"),quote = FALSE)
    i=i+1
  }
  
  
  #make ranked genelist
  {
    res.sort <- res[sort(res$log2FoldChange,decreasing = T, index.return = T)$ix,]
    res.sort <- res.sort[grep("ENSG",rownames(res.sort)),]
    
    genelist <- as.numeric(res.sort$log2FoldChange)
    names(genelist) <- as.character(lapply(strsplit(rownames(res.sort),"\\|"), function(x) x[1]))
    names(genelist) <- as.character(lapply(strsplit(names(genelist),".",fixed = TRUE), function(x) x[1]))
    
  }
  
  #make geneset from KEGG
  {
    PATH_ID_NAME <- read.csv("pathway/PATH_ID_NAME_KEGGplusHallmark.txt",header = TRUE, row.names = 1,sep = "\t",check.names=FALSE)
    #gmt <- PATH_ID_NAME[,c(4,1)] # local_modified
    gmt <- PATH_ID_NAME[,c(6,1)] # local_plusHallmarker
    colnames(gmt) <- c("ont","gene")
  }
  
  #GSEA by clusterprofile and plot
  {
    egmt2 <- GSEA(genelist,minGSSize = 0, maxGSSize = 1000,
                  pvalueCutoff = 1, TERM2GENE = gmt)
    write.table(egmt2@result,"./output/PBMC_GSEA.txt",quote = FALSE,sep = "\t")
    
    #gsea <- gseaplot2(egmt2,1)
    #pdf("Ribosome_PBMC_CRCvsNC.pdf")
    #plot(gsea)
    #dev.off()
  }
  
  #bubble plot
  {
    path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/09.Published_data/output/"
    docs <- dir(path)
    docs <- grep("_KEGG.csv",docs,value = TRUE)
    dataset_name <- as.character(lapply(strsplit(docs,"_KEGG"),function(x) x[1]))
    print(dataset_name)
    i=1
    KEGG={}
    while(i<=length(dataset_name)){
    KEGG_tmp <- read.csv(paste0("./output/",dataset_name[i],"_KEGG.csv"),header = TRUE)
    KEGG_tmp$dataset_name <- dataset_name[i]
    KEGG <- rbind(KEGG,KEGG_tmp,stringsAsFactors = F)
    i=i+1
    }
    
    forbubbleplot <- KEGG
    
    res <- KEGG[which(KEGG$dataset_name=="TCGA_ESCA_PrimaryTumor_vs_TissueNormal"),]
    filtered <- res[res$p.adjust < 0.05, ]
    top <- head(filtered[order(filtered$p.adjust,decreasing = FALSE),]$Description,21)
    top <- factor(top,levels=top)
    
    genes <- top[-2]
    genes <- c(
      "Ribosome","Platelet activation","Endocytosis","Hematopoietic cell lineage","Antigen processing and presentation",
      "Spliceosome","Cell cycle",
      "Ras signaling pathway","p53 signaling pathway","DNA replication"
    )
    
    dataset_forplot_all <- c("TCGA_COAD_PrimaryTumor_vs_TissueNormal","TCGA_STAD_PrimaryTumor_vs_TissueNormal","TCGA_BRCA_PrimaryTumor_vs_TissueNormal",
                         "TCGA_ESCA_PrimaryTumor_vs_TissueNormal","TCGA_LIHC_PrimaryTumor_vs_TissueNormal","TCGA_LUAD_PrimaryTumor_vs_TissueNormal",
                         "TCGA_LUSC_PrimaryTumor_vs_TissueNormal","TCGA_READ_PrimaryTumor_vs_TissueNormal","TCGA_THCA_PrimaryTumor_vs_TissueNormal",
                         "GSE40174_PDAC_vs_NC","GSE120663_HCC_vs_NC","GSE117623_HCC_vs_DC",
                         "exoRBase_HCC_vs_NC","exoRBase_CRC_vs_NC","exoRBase_PAAD_vs_NC","GSE133684_PDAC_vs_NC",
                         "GSE142978_HCC_vs_NC","GSE174302_CRC_vs_NC","GSE174302_ESCA_vs_NC",
                         "GSE174302_HCC_vs_NC","GSE174302_LUAD_vs_NC","GSE174302_STAD_vs_NC",
                         "GSE68086_BRCA_vs_NC","GSE68086_CRC_vs_NC","GSE68086_GBM_vs_NC",
                         "GSE68086_Hepatobiliary_vs_NC","GSE68086_NSCLC_vs_NC","GSE68086_PDAC_vs_NC","GSE89843_NSCLC_vs_NC",
                         "GSE89843_Multiple-Sclerosis_vs_NC","GSE89843_Pulmonary-Hypertension_vs_NC","GSE89843_Chronic-Pancreatitis_vs_NC",
                         "GSE89843_Epilepsy_vs_NC","GSE89843_Atherosclerosis_vs_NC","GSE89843_Angina-Pectoris_vs_NC")
    dataset_forplot <- c("TCGA_COAD_PrimaryTumor_vs_TissueNormal","TCGA_READ_PrimaryTumor_vs_TissueNormal","TCGA_STAD_PrimaryTumor_vs_TissueNormal",
                         "TCGA_BRCA_PrimaryTumor_vs_TissueNormal","TCGA_ESCA_PrimaryTumor_vs_TissueNormal","TCGA_LIHC_PrimaryTumor_vs_TissueNormal",
                         "exoRBase_HCC_vs_NC","exoRBase_CRC_vs_NC","exoRBase_PAAD_vs_NC","GSE133684_PDAC_vs_NC",
                         "GSE142978_HCC_vs_NC","GSE174302_CRC_vs_NC","GSE174302_ESCA_vs_NC",
                         "GSE174302_HCC_vs_NC","GSE174302_LUAD_vs_NC","GSE174302_STAD_vs_NC",
                         "GSE68086_BRCA_vs_NC","GSE68086_CRC_vs_NC","GSE68086_GBM_vs_NC",
                         "GSE68086_Hepatobiliary_vs_NC","GSE68086_NSCLC_vs_NC","GSE68086_PDAC_vs_NC")
    forbubbleplot_plot <- filter(forbubbleplot, Description %in% genes)
    forbubbleplot_plot <- filter(forbubbleplot_plot, dataset_name %in% dataset_forplot)
    forbubbleplot_plot[forbubbleplot_plot$p.adjust>=0.05,"p.adjust"] <- NA
    
    forbubbleplot_plot$Description <- factor(forbubbleplot_plot$Description,levels = rev(genes))
    forbubbleplot_plot$dataset_name <- factor(forbubbleplot_plot$dataset_name,levels = c("TCGA_COAD_PrimaryTumor_vs_TissueNormal","TCGA_READ_PrimaryTumor_vs_TissueNormal","TCGA_STAD_PrimaryTumor_vs_TissueNormal","TCGA_BRCA_PrimaryTumor_vs_TissueNormal",
                                                                                                  "TCGA_ESCA_PrimaryTumor_vs_TissueNormal","TCGA_LIHC_PrimaryTumor_vs_TissueNormal","TCGA_LUAD_PrimaryTumor_vs_TissueNormal",
                                                                                                  "TCGA_LUSC_PrimaryTumor_vs_TissueNormal","TCGA_THCA_PrimaryTumor_vs_TissueNormal",
                                                                                                  "GSE40174_PDAC_vs_NC","GSE120663_HCC_vs_NC","GSE117623_HCC_vs_DC",
                                                                                                  "exoRBase_HCC_vs_NC","exoRBase_CRC_vs_NC","exoRBase_PAAD_vs_NC","GSE133684_PDAC_vs_NC",
                                                                                                  "GSE142978_HCC_vs_NC","GSE174302_CRC_vs_NC","GSE174302_ESCA_vs_NC",
                                                                                                  "GSE174302_HCC_vs_NC","GSE174302_LUAD_vs_NC","GSE174302_STAD_vs_NC",
                                                                                                  "GSE68086_BRCA_vs_NC","GSE68086_CRC_vs_NC","GSE68086_GBM_vs_NC",
                                                                                                  "GSE68086_Hepatobiliary_vs_NC","GSE68086_NSCLC_vs_NC","GSE68086_PDAC_vs_NC","GSE89843_NSCLC_vs_NC",
                                                                                                  "GSE89843_Multiple-Sclerosis_vs_NC","GSE89843_Pulmonary-Hypertension_vs_NC","GSE89843_Chronic-Pancreatitis_vs_NC",
                                                                                                  "GSE89843_Epilepsy_vs_NC","GSE89843_Atherosclerosis_vs_NC","GSE89843_Angina-Pectoris_vs_NC"))
    ggplot(forbubbleplot_plot,aes(x=forbubbleplot_plot$dataset_name,y=forbubbleplot_plot$Description))+
      geom_point(aes(size=-1*log10(forbubbleplot_plot$p.adjust),color=as.numeric(forbubbleplot_plot$Count)))+
      scale_colour_gradient2(low="#162252",high="#3B5E2B",mid="#008B00",midpoint = 0)+
      geom_vline(xintercept = c(6.5,10.5,16.5),linetype = "dashed",color="black")+
      geom_vline(xintercept = 25.5,linetype = "dashed",color="grey")+
      labs(color="GeneCount",
           size=expression(-log[10](p.adjust)),
           x="dataset")+
      theme_bw()+
      theme(axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
            axis.text.x = element_text(size=rel(1.3),face="bold",colour = "black",angle = 45,vjust=1,hjust=1),
            axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
            axis.title.y = element_blank())
  }
}

#differential expression between CRC and NC in multiomics, select genes (FDR<0.05) and then enrich at pathway level by KEGG
{
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/11.variation_KEGG/")
    path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/11.variation_KEGG/"
    docs <- dir(path)
    docs <- grep("des_multiomics.csv",docs,value = TRUE, invert = TRUE)
    docs <- grep("Altpromoter",docs,value = TRUE, invert = TRUE)
    dataset_name <- docs
    print(dataset_name)
    j=1
    while(j<=length(dataset_name)){
    mat_raw <- read.table(paste0(path,dataset_name[j],"/",dataset_name[j],"_ML.txt"), header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    des <- read.csv(paste0(path,dataset_name[j],"/des_",dataset_name[j],".csv"), header = TRUE, check.names=FALSE, sep=',')
    
    samples <- des$samples
    #samples <- gsub(".","-",samples,fixed=TRUE)
    #samples <- gsub("X","",samples,fixed=TRUE)
    group <- des$group
    batch <- des$batch
    
    i=1
    mat=as.data.frame(array(dim=c(length(rownames(mat_raw)),1)))
    while (i<=length(samples)) {
      temp <- mat_raw[,which(colnames(mat_raw)==samples[i])]
      #temp <- mat_raw[,grep(samples[i],colnames(mat_raw),fixed=TRUE)]
      #print(i)
      #print(which(colnames(mat_raw)==samples[i]))
      temp <- as.data.frame(temp)
      colnames(temp) <- samples[i]
      mat <- cbind(mat,temp)
      i=i+1
    }
    
    mat <- mat[,-1]
    rownames(mat) <- rownames(mat_raw)
    
    y <- DGEList(counts=mat, samples=samples, group=group)     #规定输入格式，必须包含这几个参数
    y <- calcNormFactors(y, method="TMM")                      #归一化处理
    y <- estimateDisp(y)
    test <- exactTest(y, pair = c("negative","positive"), dispersion = "auto")
    #https://support.bioconductor.org/p/64807/
    #广义线性模型可以针对多因素的情况，在这个场景下，只比较cancer与normal用exactTest就足够了
    {
      #design <- model.matrix(~0+group)   
      #y <- estimateDisp(y, design)          #广义线性模型计算离散度（common&trended&tagwise）
      #fit <- glmFit(y, design)              #差异表达分析函数                     
      #test <- glmLRT(fit, coef=2)           #差异表达分析函数 
    }
    
    
    
    res <- topTags(test, n=nrow(mat), sort.by='none')     #输出计算差异基因的结果
    res <- cbind(res$table, baseMean=2^(res$table$logCPM)) #输出差异基因结果加上baseMean一栏
    mapped_names <- colnames(res)                   #后面都是整理+输出
    for(i in 1:ncol(res)){
      if(colnames(res)[i] == 'logFC'){
        mapped_names[i] <- 'log2FoldChange'
      }else if(colnames(res)[i] == 'PValue'){
        mapped_names[i] <- 'pvalue'
      }else if(colnames(res)[i] == 'FDR') {
        mapped_names[i] <- 'padj'
      }else{
        mapped_names[i] <- colnames(res)[i]
      }
    }
    colnames(res) <- mapped_names
    write.table(res, paste0(path,dataset_name[j],"/",dataset_name[j],"_edger_exact.csv"), sep=',', quote=FALSE, row.names=TRUE) #输出文件，规定名字
    j=j+1
  }
  }
  
  
  
  path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/11.variation_KEGG/"
  docs <- dir(path)
  docs <- grep("des_multiomics.csv",docs,value = TRUE, invert = TRUE)
  docs <- grep("Altpromoter",docs,value = TRUE, invert = TRUE)
  dataset_name <- docs
  print(dataset_name)
  i=1
  while(i<=length(dataset_name)){
    #select gene FDR < 0.05
    DE <- read.csv(paste0(path,dataset_name[i],"/",dataset_name[i],"_edger_exact.csv"),header = TRUE, row.names = 1)
    
    all_gene <- row.names(DE)
    all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
    filtered <- rownames(DE[DE$padj<0.1,])
    filtered <- as.character(lapply(strsplit(filtered,"\\."),function(x) x[1]))
    
    forenrich <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                       filters = "ensembl_gene_id",
                       values=filtered, mart= mart,useCache = FALSE)
    background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                        filters = "ensembl_gene_id",
                        values=all_gene, mart= mart,useCache = FALSE)
    
    forenrich <- na.omit(forenrich)
    background <- na.omit(background)
    
    KEGG_res <- enrichKEGG(
      forenrich$entrezgene_id,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      universe=as.character(background$entrezgene_id),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1,
      use_internal_data = FALSE
    )
    
    #dotplot(KEGG_res)
    
    write.csv(KEGG_res,paste0(path,dataset_name[i],"/",dataset_name[i],"_KEGG.csv"),quote = FALSE)
    i=i+1
  }
  
  
  #make ranked genelist
  {
    res.sort <- res[sort(res$log2FoldChange,decreasing = T, index.return = T)$ix,]
    res.sort <- res.sort[grep("ENSG",rownames(res.sort)),]
    
    genelist <- as.numeric(res.sort$log2FoldChange)
    names(genelist) <- as.character(lapply(strsplit(rownames(res.sort),"\\|"), function(x) x[1]))
    names(genelist) <- as.character(lapply(strsplit(names(genelist),".",fixed = TRUE), function(x) x[1]))
    
  }
  
  #make geneset from KEGG
  {
    PATH_ID_NAME <- read.csv("pathway/PATH_ID_NAME_KEGGplusHallmark.txt",header = TRUE, row.names = 1,sep = "\t",check.names=FALSE)
    #gmt <- PATH_ID_NAME[,c(4,1)] # local_modified
    gmt <- PATH_ID_NAME[,c(6,1)] # local_plusHallmarker
    colnames(gmt) <- c("ont","gene")
  }
  
  #GSEA by clusterprofile and plot
  {
    egmt2 <- GSEA(genelist,minGSSize = 0, maxGSSize = 1000,
                  pvalueCutoff = 1, TERM2GENE = gmt)
    write.table(egmt2@result,"./output/PBMC_GSEA.txt",quote = FALSE,sep = "\t")
    
    #gsea <- gseaplot2(egmt2,1)
    #pdf("Ribosome_PBMC_CRCvsNC.pdf")
    #plot(gsea)
    #dev.off()
  }
  
  #bubble plot
  {
    path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/11.variation_KEGG/"
    docs <- dir(path)
    docs <- grep("Altpromoter",docs,value = TRUE, invert = TRUE)
    docs <- grep("Medip",docs,value = TRUE, invert = TRUE)
    docs <- grep("Nucleosome",docs,value = TRUE, invert = TRUE)
    dataset_name <- docs
    print(dataset_name)
    i=1
    KEGG={}
    while(i<=length(dataset_name)){
      KEGG_tmp <- read.csv(paste0(path,dataset_name[i],"/",dataset_name[i],"_KEGG.csv"),header = TRUE)
      KEGG_tmp$dataset_name <- dataset_name[i]
      KEGG <- rbind(KEGG,KEGG_tmp,stringsAsFactors = F)
      i=i+1
    }
    
    forbubbleplot <- KEGG
    
    res <- KEGG[which(KEGG$dataset_name=="CNV"),]
    filtered <- res[res$p.adjust < 0.05, ]
    top <- head(filtered[order(filtered$p.adjust,decreasing = FALSE),]$Description,10)
    top <- factor(top,levels=top)
    
    genes <- top
    genes <- c(
      "Ribosome","Ribosome biogenesis in eukaryotes","Spliceosome","B cell receptor signaling pathway","Epstein-Barr virus infection",
      "Diabetic cardiomyopathy",
      "Oxidative phosphorylation","Chemical carcinogenesis - reactive oxygen species"
    )
    
    dataset_forplot_all <- c("TCGA_COAD_PrimaryTumor_vs_TissueNormal","TCGA_STAD_PrimaryTumor_vs_TissueNormal","TCGA_BRCA_PrimaryTumor_vs_TissueNormal",
                             "TCGA_ESCA_PrimaryTumor_vs_TissueNormal","TCGA_LIHC_PrimaryTumor_vs_TissueNormal","TCGA_LUAD_PrimaryTumor_vs_TissueNormal",
                             "TCGA_LUSC_PrimaryTumor_vs_TissueNormal","TCGA_READ_PrimaryTumor_vs_TissueNormal","TCGA_THCA_PrimaryTumor_vs_TissueNormal",
                             "GSE40174_PDAC_vs_NC","GSE120663_HCC_vs_NC","GSE117623_HCC_vs_DC",
                             "exoRBase_HCC_vs_NC","exoRBase_CRC_vs_NC","exoRBase_PAAD_vs_NC","GSE133684_PDAC_vs_NC",
                             "GSE142978_HCC_vs_NC","GSE174302_CRC_vs_NC","GSE174302_ESCA_vs_NC",
                             "GSE174302_HCC_vs_NC","GSE174302_LUAD_vs_NC","GSE174302_STAD_vs_NC",
                             "GSE68086_BRCA_vs_NC","GSE68086_CRC_vs_NC","GSE68086_GBM_vs_NC",
                             "GSE68086_Hepatobiliary_vs_NC","GSE68086_NSCLC_vs_NC","GSE68086_PDAC_vs_NC","GSE89843_NSCLC_vs_NC",
                             "GSE89843_Multiple-Sclerosis_vs_NC","GSE89843_Pulmonary-Hypertension_vs_NC","GSE89843_Chronic-Pancreatitis_vs_NC",
                             "GSE89843_Epilepsy_vs_NC","GSE89843_Atherosclerosis_vs_NC","GSE89843_Angina-Pectoris_vs_NC")
    dataset_forplot <- c("TCGA_COAD_PrimaryTumor_vs_TissueNormal","TCGA_READ_PrimaryTumor_vs_TissueNormal","TCGA_STAD_PrimaryTumor_vs_TissueNormal",
                         "TCGA_BRCA_PrimaryTumor_vs_TissueNormal","TCGA_ESCA_PrimaryTumor_vs_TissueNormal","TCGA_LIHC_PrimaryTumor_vs_TissueNormal",
                         "exoRBase_HCC_vs_NC","exoRBase_CRC_vs_NC","exoRBase_PAAD_vs_NC","GSE133684_PDAC_vs_NC",
                         "GSE142978_HCC_vs_NC","GSE174302_CRC_vs_NC","GSE174302_ESCA_vs_NC",
                         "GSE174302_HCC_vs_NC","GSE174302_LUAD_vs_NC","GSE174302_STAD_vs_NC",
                         "GSE68086_BRCA_vs_NC","GSE68086_CRC_vs_NC","GSE68086_GBM_vs_NC",
                         "GSE68086_Hepatobiliary_vs_NC","GSE68086_NSCLC_vs_NC","GSE68086_PDAC_vs_NC")
    forbubbleplot_plot <- filter(forbubbleplot, Description %in% genes)
    forbubbleplot_plot <- filter(forbubbleplot_plot, dataset_name %in% dataset_forplot)
    forbubbleplot_plot[forbubbleplot_plot$p.adjust>=0.05,"p.adjust"] <- NA
    
    forbubbleplot_plot$Description <- factor(forbubbleplot_plot$Description,levels = rev(genes))
    forbubbleplot_plot$dataset_name <- factor(forbubbleplot_plot$dataset_name,levels = c("TCGA_COAD_PrimaryTumor_vs_TissueNormal","TCGA_READ_PrimaryTumor_vs_TissueNormal","TCGA_STAD_PrimaryTumor_vs_TissueNormal","TCGA_BRCA_PrimaryTumor_vs_TissueNormal",
                                                                                         "TCGA_ESCA_PrimaryTumor_vs_TissueNormal","TCGA_LIHC_PrimaryTumor_vs_TissueNormal","TCGA_LUAD_PrimaryTumor_vs_TissueNormal",
                                                                                         "TCGA_LUSC_PrimaryTumor_vs_TissueNormal","TCGA_THCA_PrimaryTumor_vs_TissueNormal",
                                                                                         "GSE40174_PDAC_vs_NC","GSE120663_HCC_vs_NC","GSE117623_HCC_vs_DC",
                                                                                         "exoRBase_HCC_vs_NC","exoRBase_CRC_vs_NC","exoRBase_PAAD_vs_NC","GSE133684_PDAC_vs_NC",
                                                                                         "GSE142978_HCC_vs_NC","GSE174302_CRC_vs_NC","GSE174302_ESCA_vs_NC",
                                                                                         "GSE174302_HCC_vs_NC","GSE174302_LUAD_vs_NC","GSE174302_STAD_vs_NC",
                                                                                         "GSE68086_BRCA_vs_NC","GSE68086_CRC_vs_NC","GSE68086_GBM_vs_NC",
                                                                                         "GSE68086_Hepatobiliary_vs_NC","GSE68086_NSCLC_vs_NC","GSE68086_PDAC_vs_NC","GSE89843_NSCLC_vs_NC",
                                                                                         "GSE89843_Multiple-Sclerosis_vs_NC","GSE89843_Pulmonary-Hypertension_vs_NC","GSE89843_Chronic-Pancreatitis_vs_NC",
                                                                                         "GSE89843_Epilepsy_vs_NC","GSE89843_Atherosclerosis_vs_NC","GSE89843_Angina-Pectoris_vs_NC"))
    ggplot(forbubbleplot_plot,aes(x=forbubbleplot_plot$dataset_name,y=forbubbleplot_plot$Description))+
      geom_point(aes(size=-1*log10(forbubbleplot_plot$p.adjust),color=as.numeric(forbubbleplot_plot$Count)))+
      scale_colour_gradient2(low="#162252",high="#3B5E2B",mid="#008B00",midpoint = 0)+
      geom_vline(xintercept = c(6.5,10.5,16.5),linetype = "dashed",color="black")+
      geom_vline(xintercept = 25.5,linetype = "dashed",color="grey")+
      labs(color="GeneCount",
           size=expression(-log[10](p.adjust)),
           x="dataset")+
      theme_bw()+
      theme(axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
            axis.text.x = element_text(size=rel(1.3),face="bold",colour = "black",angle = 45,vjust=1,hjust=1),
            axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
            axis.title.y = element_blank())
    
  }
}
 
#Correlation between tissue and plasma at gene level, then enrich at pathway level by KEGG
{
  #correlation
  {
    library(pspearman)
    spearman_CI <- function(x, y, alpha = 0.05){
      rs <- cor(x, y, method = "spearman", use = "complete.obs")
      n <- sum(complete.cases(x, y))
      sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
    }
    
    output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation_at_gene_level/Plasma_vs_tumor_expression/"
    cor_methods <- "spearman"
    counts_plasma <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/plasma/TPM_intron_spanning_noMTRNA.txt",sep="\t",header = T, row.names = 1)
    counts_tissue <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/tissue/tissue_TPM_intron_spanning_noMTRNA.txt",sep="\t",header = T, row.names = 1)
    counts_PBMC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/PBMC/TPM_intron_spanning_noMTRNA.txt",sep="\t",header = T, row.names = 1)
    
    counts_plasma_tissue <- cbind(counts_plasma,counts_tissue)
    
    plasma_forcor <- counts_plasma_tissue[,c("CRC.PKU.27.pico","CRC.PKU.28.pico","CRC.PKU.29.pico","CRC.PKU.30.pico",
                                             "CRC.PKU.32.pico","CRC.PKU.34.pico","CRC.PKU.35.pico","CRC.PKU.36.pico",
                                             "CRC.PKU.37.pico","CRC.PKU.38.pico","CRC.PKU.39.pico","CRC.PKU.40.pico","CRC.PKU.41.pico")]
    colnames(plasma_forcor) <- gsub(".","-",fixed=TRUE,colnames(plasma_forcor))
    colnames(plasma_forcor) <- gsub("-pico","-plasma",fixed=TRUE,colnames(plasma_forcor))
    plasma_forcor <- t(plasma_forcor)
    tumor_forcor <- counts_plasma_tissue[,c("CRC.PKU.27.T","CRC.PKU.28.T","CRC.PKU.29.T","CRC.PKU.30.T",
                                            "CRC.PKU.32.T","CRC.PKU.34.T","CRC.PKU.35.T","CRC.PKU.36.T",
                                            "CRC.PKU.37.T","CRC.PKU.38.T","CRC.PKU.39.T","CRC.PKU.40.T","CRC.PKU.41.T")]
    colnames(tumor_forcor) <- gsub(".","-",fixed=TRUE,colnames(tumor_forcor))
    colnames(tumor_forcor) <- gsub("-T","-tumor",fixed=TRUE,colnames(tumor_forcor))
    tumor_forcor <- t(tumor_forcor)
    normal_forcor <- counts_plasma_tissue[,c("CRC.PKU.27.N","CRC.PKU.28.N","CRC.PKU.29.N","CRC.PKU.30.N",
                                             "CRC.PKU.32.N","CRC.PKU.34.N","CRC.PKU.35.N","CRC.PKU.36.N",
                                             "CRC.PKU.37.N","CRC.PKU.38.N","CRC.PKU.39.N","CRC.PKU.40.N","CRC.PKU.41.N")]
    colnames(normal_forcor) <- gsub(".","-",fixed=TRUE,colnames(normal_forcor))
    colnames(normal_forcor) <- gsub("-N","-normal",fixed=TRUE,colnames(normal_forcor))
    normal_forcor <- t(normal_forcor)
    
    #plasma vs tumor
    {
    pb <- progress_bar$new(
      format = "  Processing [:bar] :percent eta: :eta",
      total = ncol(plasma_forcor), clear = FALSE, width= 60)
    
    result_final <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result_final) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    
    result <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    
    i=1
    while(i<=(ncol(plasma_forcor))){
          r_twosided <- spearman.test(plasma_forcor[,i],tumor_forcor[,i],alternative = "two.sided", approximation = "exact")
          if (is.na(r_twosided$estimate)) {
            r <- spearman.test(plasma_forcor[,i],tumor_forcor[,i],alternative = "two.sided", approximation = "exact")
          } else if(r_twosided$estimate>0){
            r <- spearman.test(plasma_forcor[,i],tumor_forcor[,i],alternative = "greater", approximation = "exact")
          } else if (r_twosided$estimate<0) {
            r <- spearman.test(plasma_forcor[,i],tumor_forcor[,i],alternative = "less", approximation = "exact")
          } else {
            r <- spearman.test(plasma_forcor[,i],tumor_forcor[,i],alternative = "two.sided", approximation = "exact")
          }
          result_tmp <- data.frame(paste0(colnames(plasma_forcor)[i],"-plasma vs ",colnames(tumor_forcor)[i],"-tumor"),
                                   r$estimate,r$p.value,r$method,
                                   spearman_CI(plasma_forcor[,i],tumor_forcor[,i])[1],spearman_CI(plasma_forcor[,i],tumor_forcor[,i])[2],r_twosided$p.value)
        
        colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
        result <- rbind(result,result_tmp)

      if(i%%5000==0){
        result_final <- rbind(result_final,result)
        write.csv(result,paste0(output_dir,i/5000,"_Result5000.csv"))
        result <- as.data.frame(matrix(numeric(0),ncol=7))
        colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
        i=i+1
      } else if(i==(ncol(plasma_forcor))){
        result_final <- rbind(result_final,result)
        write.csv(result,paste0(output_dir,i%/%5000+1,"_Result",i%%5000,".csv"))
        i=i+1
      } else {
        i=i+1
      }
      pb$tick()
      Sys.sleep(1 / 100)
    }
    
    result_final$FDR <- p.adjust(result_final$Pvalue,"BH",n=length(result_final$Pvalue))
    
    write.csv(result_final,paste0(output_dir,"Result_final_plasma_vs_Tumor.csv"))
    
    result_final <- result_final[grep("ENSG",result_final$DataNames),]
    
    result_final$FDR <- p.adjust(result_final$Pvalue,"BH",n=length(result_final$Pvalue))
    
    write.csv(result_final,paste0(output_dir,"Result_final_plasma_vs_Tumor_ensembl.csv"))
    }
    
    #plasma vs normal
    {
    output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation_at_gene_level/Plamsa_vs_normal_expression/"
    pb <- progress_bar$new(
      format = "  Processing [:bar] :percent eta: :eta",
      total = ncol(plasma_forcor), clear = FALSE, width= 60)
    
    result_final <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result_final) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    
    result <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    
    i=1
    while(i<=(ncol(plasma_forcor))){
      r_twosided <- spearman.test(plasma_forcor[,i],normal_forcor[,i],alternative = "two.sided", approximation = "exact")
      if (is.na(r_twosided$estimate)) {
        r <- spearman.test(plasma_forcor[,i],normal_forcor[,i],alternative = "two.sided", approximation = "exact")
      } else if(r_twosided$estimate>0){
        r <- spearman.test(plasma_forcor[,i],normal_forcor[,i],alternative = "greater", approximation = "exact")
      } else if (r_twosided$estimate<0) {
        r <- spearman.test(plasma_forcor[,i],normal_forcor[,i],alternative = "less", approximation = "exact")
      } else {
        r <- spearman.test(plasma_forcor[,i],normal_forcor[,i],alternative = "two.sided", approximation = "exact")
      }
      result_tmp <- data.frame(paste0(colnames(plasma_forcor)[i],"-plasma vs ",colnames(normal_forcor)[i],"-normal"),
                               r$estimate,r$p.value,r$method,
                               spearman_CI(plasma_forcor[,i],normal_forcor[,i])[1],spearman_CI(plasma_forcor[,i],normal_forcor[,i])[2],r_twosided$p.value)
      
      colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
      result <- rbind(result,result_tmp)
      
      if(i%%5000==0){
        result_final <- rbind(result_final,result)
        write.csv(result,paste0(output_dir,i/5000,"_Result5000.csv"))
        result <- as.data.frame(matrix(numeric(0),ncol=7))
        colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
        i=i+1
      } else if(i==(ncol(plasma_forcor))){
        result_final <- rbind(result_final,result)
        write.csv(result,paste0(output_dir,i%/%5000+1,"_Result",i%%5000,".csv"))
        i=i+1
      } else {
        i=i+1
      }
      pb$tick()
      Sys.sleep(1 / 100)
    }
    
    result_final$FDR <- p.adjust(result_final$Pvalue,"BH",n=length(result_final$Pvalue))
    
    write.csv(result_final,paste0(output_dir,"Result_final_plasma_vs_Normal.csv"))
    
    result_final <- result_final[grep("ENSG",result_final$DataNames),]
    
    result_final$FDR <- p.adjust(result_final$Pvalue,"BH",n=length(result_final$Pvalue))
    
    write.csv(result_final,paste0(output_dir,"Result_final_plasma_vs_Normal_ensembl.csv"))
    
    }
    
    plasma_PBMC_CRC_forcor <- counts_plasma[,c("CRC.PKU.27.pico","CRC.PKU.28.pico","CRC.PKU.29.pico",
                                             "CRC.PKU.30.pico","CRC.PKU.31.pico","CRC.PKU.32.pico")]
    colnames(plasma_PBMC_CRC_forcor) <- gsub(".","-",fixed=TRUE,colnames(plasma_PBMC_CRC_forcor))
    colnames(plasma_PBMC_CRC_forcor) <- gsub("-pico","-plasma",fixed=TRUE,colnames(plasma_PBMC_CRC_forcor))
    plasma_PBMC_CRC_forcor <- t(plasma_PBMC_CRC_forcor)
    
    plasma_PBMC_NC_forcor <- counts_plasma[,c("NC.PKU.10.pico","NC.PKU.11.pico","NC.PKU.12.pico",
                                              "NC.PKU.13.pico","NC.PKU.14.pico","NC.PKU.9.pico")]
    colnames(plasma_PBMC_NC_forcor) <- gsub(".","-",fixed=TRUE,colnames(plasma_PBMC_NC_forcor))
    colnames(plasma_PBMC_NC_forcor) <- gsub("-pico","-plasma",fixed=TRUE,colnames(plasma_PBMC_NC_forcor))
    plasma_PBMC_NC_forcor <- t(plasma_PBMC_NC_forcor)
    
    PBMC_CRC_forcor <- counts_PBMC[,c("CRC.PKU.27.PBMC","CRC.PKU.28.PBMC","CRC.PKU.29.PBMC",
                                      "CRC.PKU.30.PBMC","CRC.PKU.31.PBMC","CRC.PKU.32.PBMC")]
    colnames(PBMC_CRC_forcor) <- gsub(".","-",fixed=TRUE,colnames(PBMC_CRC_forcor))
    PBMC_CRC_forcor <- t(PBMC_CRC_forcor)
    
    PBMC_NC_forcor <- counts_PBMC[,c("NC.PKU.10.PBMC","NC.PKU.11.PBMC","NC.PKU.12.PBMC",
                                     "NC.PKU.13.PBMC","NC.PKU.14.PBMC","NC.PKU.9.PBMC")]
    colnames(PBMC_NC_forcor) <- gsub(".","-",fixed=TRUE,colnames(PBMC_NC_forcor))
    PBMC_NC_forcor <- t(PBMC_NC_forcor)
    
    #plasma vs PBMC CRC
    {
      output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation_at_gene_level/Plamsa_vs_PBMC_CRC_expression/"
      pb <- progress_bar$new(
        format = "  Processing [:bar] :percent eta: :eta",
        total = ncol(plasma_PBMC_CRC_forcor), clear = FALSE, width= 60)
      
      result_final <- as.data.frame(matrix(numeric(0),ncol=7))
      colnames(result_final) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
      
      result <- as.data.frame(matrix(numeric(0),ncol=7))
      colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
      
      i=1
      while(i<=(ncol(plasma_PBMC_CRC_forcor))){
        r_twosided <- spearman.test(plasma_PBMC_CRC_forcor[,i],PBMC_CRC_forcor[,i],alternative = "two.sided", approximation = "exact")
        if (is.na(r_twosided$estimate)) {
          r <- spearman.test(plasma_PBMC_CRC_forcor[,i],PBMC_CRC_forcor[,i],alternative = "two.sided", approximation = "exact")
        } else if(r_twosided$estimate>0){
          r <- spearman.test(plasma_PBMC_CRC_forcor[,i],PBMC_CRC_forcor[,i],alternative = "greater", approximation = "exact")
        } else if (r_twosided$estimate<0) {
          r <- spearman.test(plasma_PBMC_CRC_forcor[,i],PBMC_CRC_forcor[,i],alternative = "less", approximation = "exact")
        } else {
          r <- spearman.test(plasma_PBMC_CRC_forcor[,i],PBMC_CRC_forcor[,i],alternative = "two.sided", approximation = "exact")
        }
        result_tmp <- data.frame(paste0(colnames(plasma_PBMC_CRC_forcor)[i],"-plasma vs ",colnames(PBMC_CRC_forcor)[i],"-normal"),
                                 r$estimate,r$p.value,r$method,
                                 spearman_CI(plasma_PBMC_CRC_forcor[,i],PBMC_CRC_forcor[,i])[1],spearman_CI(plasma_PBMC_CRC_forcor[,i],PBMC_CRC_forcor[,i])[2],r_twosided$p.value)
        
        colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
        result <- rbind(result,result_tmp)
        
        if(i%%5000==0){
          result_final <- rbind(result_final,result)
          write.csv(result,paste0(output_dir,i/5000,"_Result5000.csv"))
          result <- as.data.frame(matrix(numeric(0),ncol=7))
          colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
          i=i+1
        } else if(i==(ncol(plasma_PBMC_CRC_forcor))){
          result_final <- rbind(result_final,result)
          write.csv(result,paste0(output_dir,i%/%5000+1,"_Result",i%%5000,".csv"))
          i=i+1
        } else {
          i=i+1
        }
        pb$tick()
        Sys.sleep(1 / 100)
      }
      
      result_final$FDR <- p.adjust(result_final$Pvalue,"BH",n=length(result_final$Pvalue))
      
      write.csv(result_final,paste0(output_dir,"Result_final_plasma_vs_PBMC_CRC.csv"))
      
      result_final <- result_final[grep("ENSG",result_final$DataNames),]
      
      result_final$FDR <- p.adjust(result_final$Pvalue,"BH",n=length(result_final$Pvalue))
      
      write.csv(result_final,paste0(output_dir,"Result_final_plasma_vs_PBMC_CRC_ensembl.csv"))
      
    }
    
    #plasma vs PBMC NC
    {
      output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation_at_gene_level/Plamsa_vs_PBMC_NC_expression/"
      pb <- progress_bar$new(
        format = "  Processing [:bar] :percent eta: :eta",
        total = ncol(plasma_PBMC_NC_forcor), clear = FALSE, width= 60)
      
      result_final <- as.data.frame(matrix(numeric(0),ncol=7))
      colnames(result_final) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
      
      result <- as.data.frame(matrix(numeric(0),ncol=7))
      colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
      
      i=1
      while(i<=(ncol(plasma_PBMC_NC_forcor))){
        r_twosided <- spearman.test(plasma_PBMC_NC_forcor[,i],PBMC_NC_forcor[,i],alternative = "two.sided", approximation = "exact")
        if (is.na(r_twosided$estimate)) {
          r <- spearman.test(plasma_PBMC_NC_forcor[,i],PBMC_NC_forcor[,i],alternative = "two.sided", approximation = "exact")
        } else if(r_twosided$estimate>0){
          r <- spearman.test(plasma_PBMC_NC_forcor[,i],PBMC_NC_forcor[,i],alternative = "greater", approximation = "exact")
        } else if (r_twosided$estimate<0) {
          r <- spearman.test(plasma_PBMC_NC_forcor[,i],PBMC_NC_forcor[,i],alternative = "less", approximation = "exact")
        } else {
          r <- spearman.test(plasma_PBMC_NC_forcor[,i],PBMC_NC_forcor[,i],alternative = "two.sided", approximation = "exact")
        }
        result_tmp <- data.frame(paste0(colnames(plasma_PBMC_NC_forcor)[i],"-plasma vs ",colnames(PBMC_NC_forcor)[i],"-normal"),
                                 r$estimate,r$p.value,r$method,
                                 spearman_CI(plasma_PBMC_NC_forcor[,i],PBMC_NC_forcor[,i])[1],spearman_CI(plasma_PBMC_NC_forcor[,i],PBMC_NC_forcor[,i])[2],r_twosided$p.value)
        
        colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
        result <- rbind(result,result_tmp)
        
        if(i%%5000==0){
          result_final <- rbind(result_final,result)
          write.csv(result,paste0(output_dir,i/5000,"_Result5000.csv"))
          result <- as.data.frame(matrix(numeric(0),ncol=7))
          colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
          i=i+1
        } else if(i==(ncol(plasma_PBMC_NC_forcor))){
          result_final <- rbind(result_final,result)
          write.csv(result,paste0(output_dir,i%/%5000+1,"_Result",i%%5000,".csv"))
          i=i+1
        } else {
          i=i+1
        }
        pb$tick()
        Sys.sleep(1 / 100)
      }
      
      result_final$FDR <- p.adjust(result_final$Pvalue,"BH",n=length(result_final$Pvalue))
      
      write.csv(result_final,paste0(output_dir,"Result_final_plasma_vs_PBMC_NC.csv"))
      
      result_final <- result_final[grep("ENSG",result_final$DataNames),]
      
      result_final$FDR <- p.adjust(result_final$Pvalue,"BH",n=length(result_final$Pvalue))
      
      write.csv(result_final,paste0(output_dir,"Result_final_plasma_vs_PBMC_NC_ensembl.csv"))
      
    }
    
  }
  
  #KEGG
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation_at_gene_level")
    specimen_correlation_Pvalue_0.05 <- read.csv("20210804_specimen_correlation.csv",header = TRUE)
    
    tumor_all <- specimen_correlation_Pvalue_0.05[which(specimen_correlation_Pvalue_0.05$Type=="Plasma vs tumor tissue"),]
    tumor_pos <- tumor_all[tumor_all$R>0,]
    tumor_neg <- tumor_all[tumor_all$R<0,]
    
    normal_all <- specimen_correlation_Pvalue_0.05[which(specimen_correlation_Pvalue_0.05$Type=="Plasma vs normal tissue"),]
    normal_pos <- normal_all[normal_all$R>0,]
    normal_neg <- normal_all[normal_all$R<0,]
    
    Patients_PBMC_all <- specimen_correlation_Pvalue_0.05[which(specimen_correlation_Pvalue_0.05$Type=="Patients plasma vs PBMC"),]
    Patients_PBMC_pos <- Patients_PBMC_all[Patients_PBMC_all$R>0,]
    Patients_PBMC_neg <- Patients_PBMC_all[Patients_PBMC_all$R<0,]
    
    Healthy_PBMC_all <- specimen_correlation_Pvalue_0.05[which(specimen_correlation_Pvalue_0.05$Type=="Healthy plasma vs PBMC"),]
    Healthy_PBMC_pos <- Healthy_PBMC_all[Healthy_PBMC_all$R>0,]
    Healthy_PBMC_neg <- Healthy_PBMC_all[Healthy_PBMC_all$R<0,]
    
    correlation_summary <- c("postive correlated genes between plasma and tumor"=nrow(tumor_pos),
                             "negative correlated genes between plasma and tumor"=nrow(tumor_neg),
                             "postive correlated genes between plasma and normal"=nrow(normal_pos),
                             "negative correlated genes between plasma and normal"=nrow(normal_neg),
                             "postive correlated genes between patients plasma and PBMC"=nrow(Patients_PBMC_pos),
                             "negative correlated genes between patients plasma and PBMC"=nrow(Patients_PBMC_neg),
                             "postive correlated genes between healthys plasma and PBMC"=nrow(Healthy_PBMC_pos),
                             "negative correlated genes between healthys plasma and PBMC"=nrow(Healthy_PBMC_neg))
    
    correlation_barplot <- as.data.frame(correlation_summary)
    correlation_barplot$name <- rownames(correlation_barplot)
    correlation_barplot$name <- factor(correlation_barplot$name,levels=names(correlation_summary))
    ggplot(correlation_barplot,aes(x=name,y=correlation_summary,fill=name))+
      geom_bar(stat = "identity",colour = "black")+
      geom_text(aes(x=name,y=correlation_summary+40,label=correlation_summary),size = 4,angle = 0)+
      scale_y_continuous(expand = c(0,0),limits = c(0,2000))+
      scale_fill_manual(values = c("black","grey","black","grey","black","grey","black","grey"))+
      geom_vline(xintercept = c(2.5,4.5,6.5), linetype="dashed",lwd = 0.5, color = "grey")+
      xlab("")+
      ylab("Number")+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(15,5,10,5),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.y = element_line(color = "black"),
        #axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))
    
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    
    ensembl_id <- tumor_pos$gene_id
    ensembl_id <- tumor_neg$gene_id
    ensembl_id <- normal_pos$gene_id
    ensembl_id <- normal_neg$gene_id
    ensembl_id <- Patients_PBMC_pos$gene_id
    ensembl_id <- Patients_PBMC_neg$gene_id
    ensembl_id <- Healthy_PBMC_pos$gene_id
    ensembl_id <- Healthy_PBMC_neg$gene_id
    {
    forenrich <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                            filters = "ensembl_gene_id",
                            values=ensembl_id, mart= mart,useCache = FALSE)
    forenrich <- na.omit(forenrich)
      
    KEGG_res <- enrichKEGG(
      forenrich$entrezgene_id,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 0.1,
      pAdjustMethod = "BH",
      universe,
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
    
    dotplot(KEGG_res)
    }
  }
}

#select highly correlated genes between plasma and tumor, then do wilcox test and boxplot
{
  ##prepare gene expression matrix
  #plasma(paired with tissue)
  plasma_forcor <- counts_plasma_tissue[,c("CRC.PKU.27.pico","CRC.PKU.28.pico","CRC.PKU.29.pico","CRC.PKU.30.pico",
                                             "CRC.PKU.32.pico","CRC.PKU.34.pico","CRC.PKU.35.pico","CRC.PKU.36.pico",
                                             "CRC.PKU.37.pico","CRC.PKU.38.pico","CRC.PKU.39.pico","CRC.PKU.40.pico","CRC.PKU.41.pico")]
  colnames(plasma_forcor) <- gsub(".","-",fixed=TRUE,colnames(plasma_forcor))
  colnames(plasma_forcor) <- gsub("-pico","-plasma",fixed=TRUE,colnames(plasma_forcor))
  #tumor and normal tissue (paired with plasma)
  tumor_forcor <- counts_plasma_tissue[,c("CRC.PKU.27.T","CRC.PKU.28.T","CRC.PKU.29.T","CRC.PKU.30.T",
                                          "CRC.PKU.32.T","CRC.PKU.34.T","CRC.PKU.35.T","CRC.PKU.36.T",
                                          "CRC.PKU.37.T","CRC.PKU.38.T","CRC.PKU.39.T","CRC.PKU.40.T","CRC.PKU.41.T")]
  colnames(tumor_forcor) <- gsub(".","-",fixed=TRUE,colnames(tumor_forcor))
  colnames(tumor_forcor) <- gsub("-T","-tumor",fixed=TRUE,colnames(tumor_forcor))

  normal_forcor <- counts_plasma_tissue[,c("CRC.PKU.27.N","CRC.PKU.28.N","CRC.PKU.29.N","CRC.PKU.30.N",
                                           "CRC.PKU.32.N","CRC.PKU.34.N","CRC.PKU.35.N","CRC.PKU.36.N",
                                           "CRC.PKU.37.N","CRC.PKU.38.N","CRC.PKU.39.N","CRC.PKU.40.N","CRC.PKU.41.N")]
  colnames(normal_forcor) <- gsub(".","-",fixed=TRUE,colnames(normal_forcor))
  colnames(normal_forcor) <- gsub("-N","-normal",fixed=TRUE,colnames(normal_forcor))
  
  ##prepare gene expression matrix
  #plasma(paired with PBMC)
  plasma_PBMC_CRC_forcor <- counts_plasma[,c("CRC.PKU.27.pico","CRC.PKU.28.pico","CRC.PKU.29.pico",
                                             "CRC.PKU.30.pico","CRC.PKU.31.pico","CRC.PKU.32.pico")]
  colnames(plasma_PBMC_CRC_forcor) <- gsub(".","-",fixed=TRUE,colnames(plasma_PBMC_CRC_forcor))
  colnames(plasma_PBMC_CRC_forcor) <- gsub("-pico","-plasma",fixed=TRUE,colnames(plasma_PBMC_CRC_forcor))

  plasma_PBMC_NC_forcor <- counts_plasma[,c("NC.PKU.10.pico","NC.PKU.11.pico","NC.PKU.12.pico",
                                            "NC.PKU.13.pico","NC.PKU.14.pico","NC.PKU.9.pico")]
  colnames(plasma_PBMC_NC_forcor) <- gsub(".","-",fixed=TRUE,colnames(plasma_PBMC_NC_forcor))
  colnames(plasma_PBMC_NC_forcor) <- gsub("-pico","-plasma",fixed=TRUE,colnames(plasma_PBMC_NC_forcor))

  #PBMC (paired with plasma)
  PBMC_CRC_forcor <- counts_PBMC[,c("CRC.PKU.27.PBMC","CRC.PKU.28.PBMC","CRC.PKU.29.PBMC",
                                    "CRC.PKU.30.PBMC","CRC.PKU.31.PBMC","CRC.PKU.32.PBMC")]
  colnames(PBMC_CRC_forcor) <- gsub(".","-",fixed=TRUE,colnames(PBMC_CRC_forcor))

  
  PBMC_NC_forcor <- counts_PBMC[,c("NC.PKU.10.PBMC","NC.PKU.11.PBMC","NC.PKU.12.PBMC",
                                   "NC.PKU.13.PBMC","NC.PKU.14.PBMC","NC.PKU.9.PBMC")]
  colnames(PBMC_NC_forcor) <- gsub(".","-",fixed=TRUE,colnames(PBMC_NC_forcor))

  
  {
  #input gene ids
  gene_id <- "ENSG00000111796"
  #outdir
  outdir <- paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation_at_gene_level/genes/",gene_id)
  #tissue correlation plot
  plasma <- plasma_forcor[grep(gene_id,rownames(plasma_forcor)),]
  colnames(plasma) <- gsub("-plasma","",fixed=TRUE,colnames(plasma))
  rownames(plasma) <- "plasma"
  tumor <- tumor_forcor[grep(gene_id,rownames(tumor_forcor)),]
  colnames(tumor) <- gsub("-tumor","",fixed=TRUE,colnames(tumor))
  rownames(tumor) <- "tumor"
  normal <- normal_forcor[grep(gene_id,rownames(normal_forcor)),]
  colnames(normal) <- gsub("-normal","",fixed=TRUE,colnames(normal))
  rownames(normal) <- "normal"
  
  forcor <- as.data.frame(t(rbind(plasma,tumor,normal)))
  
  Plasma_vs_Tumor_tissue <- ggscatter(forcor, x = "tumor", y = "plasma", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      plot.margin = unit(c(10,10,10,10),"pt"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
  
  Plasma_vs_Normal_tissue <- ggscatter(forcor, x = "normal", y = "plasma", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      plot.margin = unit(c(10,10,10,10),"pt"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
  
  #PBMC correlation plot
  plasma_PBMC_CRC <- plasma_PBMC_CRC_forcor[grep(gene_id,rownames(plasma_PBMC_CRC_forcor)),]
  colnames(plasma_PBMC_CRC) <- gsub("-plasma","",fixed=TRUE,colnames(plasma_PBMC_CRC))
  rownames(plasma_PBMC_CRC) <- "plasma_PBMC_CRC"
  
  plasma_PBMC_NC <- plasma_PBMC_NC_forcor[grep(gene_id,rownames(plasma_PBMC_NC_forcor)),]
  colnames(plasma_PBMC_NC) <- gsub("-plasma","",fixed=TRUE,colnames(plasma_PBMC_NC))
  rownames(plasma_PBMC_NC) <- "plasma_PBMC_NC"
  
  PBMC_CRC <- PBMC_CRC_forcor[grep(gene_id,rownames(PBMC_CRC_forcor)),]
  colnames(PBMC_CRC) <- gsub("-PBMC","",fixed=TRUE,colnames(PBMC_CRC))
  rownames(PBMC_CRC) <- "PBMC_CRC"
  
  PBMC_NC <- PBMC_NC_forcor[grep(gene_id,rownames(PBMC_NC_forcor)),]
  colnames(PBMC_NC) <- gsub("-PBMC","",fixed=TRUE,colnames(PBMC_NC))
  rownames(PBMC_NC) <- "PBMC_NC"
  
  forcor_PBMC_CRC <- as.data.frame(t(rbind(plasma_PBMC_CRC,PBMC_CRC)))
  forcor_PBMC_NC <- as.data.frame(t(rbind(plasma_PBMC_NC,PBMC_NC)))
  
  Plasma_vs_Tumor_PBMC <- ggscatter(forcor_PBMC_CRC, x = "PBMC_CRC", y = "plasma_PBMC_CRC", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      plot.margin = unit(c(10,10,10,10),"pt"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
  
  Plasma_vs_Normal_PBMC <- ggscatter(forcor_PBMC_NC, x = "PBMC_NC", y = "plasma_PBMC_NC", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 6))+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      plot.margin = unit(c(10,10,10,10),"pt"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
  
  ggsave(Plasma_vs_Tumor_tissue,filename = paste0(gene_id,"_Plasma_vs_Tumor_tissue",".pdf"),path = outdir)
  ggsave(Plasma_vs_Normal_tissue,filename = paste0(gene_id,"_Plasma_vs_Normal_tissue",".pdf"),path = outdir)
  ggsave(Plasma_vs_Tumor_PBMC,filename = paste0(gene_id,"_Plasma_vs_Tumor_PBMC",".pdf"),path = outdir)
  ggsave(Plasma_vs_Normal_PBMC,filename = paste0(gene_id,"_Plasma_vs_Normal_PBMC",".pdf"),path = outdir)
  
  boxplot_plasma_tissue <- melt(forcor)
  boxplot_plasma_tissue[,1] <- gsub("plasma","plasma_patient",boxplot_plasma_tissue[,1])
  boxplot_plasma_tissue$sample <- rep(rownames(forcor),ncol(forcor))
  
  boxplot_plasma_PBMC_NC <- melt(forcor_PBMC_NC)
  boxplot_plasma_PBMC_NC[,1] <- gsub("plasma_PBMC_NC","plasma_healthy",boxplot_plasma_PBMC_NC[,1])
  boxplot_plasma_PBMC_NC$sample <- rep(rownames(forcor_PBMC_NC),ncol(forcor_PBMC_NC))
  
  boxplot_plasma_PBMC_CRC <- melt(forcor_PBMC_CRC)
  boxplot_plasma_PBMC_CRC[,1] <- gsub("plasma_PBMC_CRC","plasma_patient",boxplot_plasma_PBMC_CRC[,1])
  boxplot_plasma_PBMC_CRC$sample <- rep(rownames(forcor_PBMC_CRC),ncol(forcor_PBMC_CRC))
  boxplot_plasma_PBMC_CRC <- boxplot_plasma_PBMC_CRC[c(5,7:12),]
  
  forboxplot <- rbind(boxplot_plasma_tissue,boxplot_plasma_PBMC_NC,boxplot_plasma_PBMC_CRC)
  
  colnames(forboxplot) <- c("Type","TPM","Sample")
  
  forboxplot$Type <- factor(forboxplot$Type,levels=c("plasma_healthy","plasma_patient","PBMC_NC","PBMC_CRC","normal","tumor"))
  p <- ggplot(forboxplot,aes(x=Type,y=log2(TPM+1),fill=Type))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
    geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
    scale_fill_brewer(palette="Blues") +
    geom_vline(xintercept=c(2.5,4.5),linetype="dashed",lwd = 1, color = "grey")+
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
  
  trend <- mean(subset(forboxplot,Type=="normal")$TPM)-mean(subset(forboxplot,Type=="tumor")$TPM)
  
  if(trend>0){
    p <- p+stat_compare_means(comparisons = list(c("plasma_healthy","plasma_patient"),c("PBMC_NC","PBMC_CRC"),c("normal","tumor")),
                              method = "wilcox.test",
                              method.args = list(alternative = "greater",exact = TRUE),
                              label = "p.signif"
    )+labs(x="",y="log2(TPM+1)",title=paste0(gene_id,"\nwilcox.test.greater"), face="bold",fill="Batch",width=10,height=8)
    ggsave(p,filename = paste0(gene_id,"_greater.pdf"),path = outdir)
  } else {
    p <- p+stat_compare_means(comparisons = list(c("plasma_healthy","plasma_patient"),c("PBMC_NC","PBMC_CRC"),c("normal","tumor")),
                              method = "wilcox.test",
                              method.args = list(alternative = "less",exact = TRUE),
                              label = "p.signif"
    )+labs(x="",y="log2(TPM+1)",title=paste0(gene_id,"\nwilcox.test.less"), face="bold",fill="Batch",width=10,height=8)
    ggsave(p,filename = paste0(gene_id,"_less.pdf"),path=outdir)
  }
  }
  #
}

##deconvolution by cibersortx
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
    library(biomaRt) ## ID convert
    library(org.Hs.eg.db) ## annotation for human
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/GTEx/second/")
    signature <- read.csv("signature_matrix_DE_top20_TPM_tissue_median_forplot.txt",sep="\t",header = FALSE,row.names = 1)
    mixture <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/tissue/tissue_TPM_intron_spanning_noMTRNA.txt",sep="\t",header = T, row.names = 1)
    
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
    write.csv(pathway_gene_count,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/tissue/tissue_forGTEx_signature_second.csv")
    #write.csv(signature,"signature_matrix_DE_top50_TPM_rmdup_ENSG.csv",quote = F)
  }
  
  #Cibersort_locally
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/GTEx/forcibersortx/")
    source('Cibersort.R')  
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
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/tissue")
    composition <- read.csv("CIBERSORTx_Job14_Results.csv",header = TRUE)
    
    library(reshape)
    stackplot <- melt(composition[,1:8])
    
    #boxplot
    {
      stackplot$variable <- gsub("\\.","-", stackplot$variable)
      stackplot$Type <- as.character(lapply(strsplit(as.character(stackplot$Mixture),".",fixed = TRUE),function(x) x[3])) 
      stackplot$Type <- gsub("T","Tumor Tissue", stackplot$Type)
      stackplot$Type <- gsub("N","Normal Tissue", stackplot$Type)
      
      my_comparisons <- list(c("Normal Tissue","Tumor Tissue"))
      p <- ggplot(stackplot[which(stackplot$variable=="Blood-Vessel"),],aes(x=Type,y=value,fill=Type))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
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
    comp_tmp <- stackplot[which(stackplot$variable=="Blood-Vessel"),]
    
    stackplot$Average <- factor(stackplot$variable,levels = unique(sort(stackplot$Average)))
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
    mixture <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/PBMC/TPM_intron_spanning_noMTRNA.txt",sep="\t",header = T, row.names = 1)
    signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/LM22/LM22.txt",sep="\t",header = T, row.names = 1)
    
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
        rownames(temp) <- gene_symbol
        colnames(temp) <- colnames(mixture)
      } else {
        #temp <- mixture[which(rownames(mixture)==target),]
        #temp <- mixture[grep(target,rownames(mixture),fixed=TRUE),]  #for ensg
        temp <- mixture[grep(target,rownames(mixture),fixed=TRUE),] #for gene symbol
        rownames(temp) <- gene_symbol
      }
      pathway_gene_count <- rbind(pathway_gene_count,temp)
      j=j+1
    }
    
    pathway_gene_count$symbol <- annotations$hgnc_symbol
    write.csv(pathway_gene_count,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/PBMC/PBMC_forLM22.csv")
  }
}


#Figure2 expression pattern

#wilcox
{
#### ribosomal protein mRNA selection workflow
library(clusterProfiler) ## enrichment
library(biomaRt) ## ID convert
library(org.Hs.eg.db) ## annotation for human
library(enrichplot) ##gseaplot
#pathway gene accession
{
  hsa_kegg <- clusterProfiler::download_KEGG("hsa")
  names(hsa_kegg)
  PATH2NAME <- hsa_kegg$KEGGPATHID2NAME
  PATH2ID <- hsa_kegg$KEGGPATHID2EXTID
  PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by="from")
  colnames(PATH_ID_NAME) <- c("KEGGID", "ENTREZID", "DESCRPTION")
  #ID conversion
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  entrezgene <- PATH_ID_NAME$ENTREZID
  ensembl_gene_id<- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                          filters = "entrezgene_id",
                          values=entrezgene, mart= mart,useCache = FALSE)
  PATH_ID_NAME <- merge(PATH_ID_NAME, ensembl_gene_id, by.x= "ENTREZID",by.y= "entrezgene_id")
  Ribosome <- subset(PATH_ID_NAME,DESCRPTION=="Ribosome")
  Ribosome_biogensis <- subset(PATH_ID_NAME,DESCRPTION=="Ribosome biogenesis in eukaryotes")
  Sulfur_metabolism <- subset(PATH_ID_NAME,DESCRPTION=="Sulfur metabolism")
  Malaria <- subset(PATH_ID_NAME,DESCRPTION=="Malaria")
  kappa <- subset(PATH_ID_NAME,DESCRPTION=="NF-kappa B signaling pathway")
  IL17 <- subset(PATH_ID_NAME,DESCRPTION=="IL-17 signaling pathway")
}

## boxpot for one pathway(need to manually specify "less","greater" or "two.sided")
#all comparison take healthy control as positive, which means "greater" is greater in healthy, "less" is less in healthy. 
{
  test_direction <- "less"
  comparison <- "NCvsCRC_STAD"
  { 
    setwd(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/",comparison,"/",test_direction))
    counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/multiomics_TPM_filtered.txt",sep="\t",header = T, row.names = 1)
    counts = log2(counts+1)
    positive <- c("NC")
    negative <- c("CRC","STAD")
    type <- c("NC","CRC","STAD") # for plot
    my_comparisons <- list(c("NC","CRC"),c("NC","STAD")) # for plot
    
    candidate <- grep("mRNA|lncRNA",rownames(counts),value = TRUE)
    candidate <- unique(candidate)

    n=1
    while(n<=length(negative)){
      
      # set progress
      pb <- progress_bar$new(
        format = "  Processing [:bar] :percent eta: :eta",
        total = length(candidate), clear = FALSE, width= 60) 
      
      j=1
      test={}
      type <- c(positive[1],negative[n])
      while(j<=length(candidate)){
        target <- candidate[j]
        if(length(grep(target,rownames(counts)))==0) {
          j=j+1
        } else {
          i=1
          gene={}
          while(i<=length(type)){
            temp <- counts[grep(target,rownames(counts),fixed=TRUE)[1],grep(type[i],colnames(counts),fixed=TRUE)]
            temp.t <- t(temp)
            colnames(temp.t) <- c("counts")
            temp.df <- as.data.frame(temp.t)
            temp.df$CancerType<-rep(type[i],times=length(temp.df$counts))
            temp.df$gene<-rep(target,times=length(temp.df$counts))
            gene <- rbind(gene,temp.df)
            i=i+1
          }
          trend <- mean(subset(gene,CancerType==positive[1])$counts)-mean(subset(gene,CancerType==negative[n])$counts)
          #wilcox.test       
          if(trend>0){
            t <- wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="greater" )
          } else if(trend<0){
            t <- wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="less" )
          } else {
            t <- wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="two.sided" )
          }
          t$gene <- rep(target,times=length(t$p.value)) 
          t$log2FC <- rep(trend,times=length(t$p.value)) 
          t$positive_vs_negative <- rep(paste0(positive[1],"_vs_",negative[n]),times=length(t$p.value))
          test <- rbind(test,t)
          j=j+1
        }
        pb$tick()
        Sys.sleep(1 / 100)
      }
      t.df={}
      test <- as.data.frame(test)
      t.df$gene <- unlist(test$gene)
      t.df$positive_vs_negative <- unlist(test$positive_vs_negative)
      t.df$p.value <- unlist(test$p.value)
      t.df$alternative <- unlist(test$alternative)
      t.df$method <- unlist(test$method)
      t.df$log2FC <- unlist(test$log2FC)
      t.df <- as.data.frame(t.df) %>% filter(p.value < 0.05,alternative==test_direction)
      write.csv(t.df,paste0(positive[1],"_vs_",negative[n],".csv"))   #输出csv文件
      y <- t.df$gene
      if(n==1){
        x <- t.df$gene
      } else {
        x <- intersect(x,y)
        print(length(x))
      }
      n=n+1
    }
    #x即为positive相对于所有negative具有秩和检验显著性的特征
    print(x)
    write.csv(x,paste0(positive[1],"_specific_feature.csv"))
  }
  
  #boxplot
  candidate <- x
  {
    j=1
    {
      test={}
      while(j<=length(as.matrix(candidate))){
        target <- as.matrix(candidate)[j]
        type <- type #keep consistant with wilcox select step
        i=1
        gene={}
        while(i<=length(type)){
          temp <- counts[grep(target,rownames(counts),fixed=TRUE)[1],grep(type[i],colnames(counts),fixed=TRUE)]
          temp.t <- t(temp)
          colnames(temp.t) <- c("counts")
          temp.df <- as.data.frame(temp.t)
          temp.df$CancerType<-rep(type[i],times=length(temp.df$counts))
          temp.df$microbe<-rep(target,times=length(temp.df$counts))
          gene <- rbind(gene,temp.df)
          i=i+1
        }
        
        my_comparisons <- my_comparisons #keep consistant with selection step
        gene$CancerType <- factor(gene$CancerType,levels=type)
        p <- ggplot(gene,aes(x=CancerType,y=counts,fill=CancerType))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
          geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
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
        
        m=test_direction
        if(m=="greater"){
          p <- p+stat_compare_means(comparisons = my_comparisons,
                                    method = "wilcox.test",
                                    method.args = list(alternative = "greater"),
                                    label = "p.signif"
          )+labs(x="",y="log2(TPM+1)",title=paste0(target,"\nwilcox.test.greater"), face="bold",fill="Type",width=10,height=8)
          ggsave(p,filename = paste0(target,"_greater.pdf"),path = "./")
        } else if(m=="less") {
          p <- p+stat_compare_means(comparisons = my_comparisons,
                                    method = "wilcox.test",
                                    method.args = list(alternative = "less"),
                                    label = "p.signif"
          )+labs(x="",y="log2(TPM+1)",title=paste0(target,"\nwilcox.test.less"), face="bold",fill="Type",width=10,height=8)
          ggsave(p,filename = paste0(target,"_less.pdf"),path="./")
        } else if(m=="two.sided") {
          p <- p+stat_compare_means(comparisons = my_comparisons,
                                    method = "wilcox.test",
                                    method.args = list(alternative = "two.sided"),
                                    label = "p.signif"
          )+labs(x="",y="log2(TPM+1)",title=paste0(target,"\nwilcox.test.two.sided"), face="bold",fill="Type",width=10,height=8)
          ggsave(p,filename = paste0(target,"_less.pdf"),path="./")
        } else {
          p <- p+stat_compare_means(comparisons = my_comparisons,
                                    method = "wilcox.test",
                                    method.args = list(alternative = "two.sided"),
                                    label = "p.signif"
          )+labs(x="",y="log2(TPM+1)",title=paste0(target,"\nwilcox.test.two.sided"), face="bold",fill="Type",width=10,height=8)
          ggsave(p,filename = paste0(target,"_less.pdf"),path="./")
        }
        j=j+1
      }
    }
  }
  
}
}

#wilcox
#greater
{
  CRCvsNC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/CRCvsSTAD_NC/greater/CRC_vs_NC.csv",header = TRUE,row.names = 1)
  CRCvsSTAD <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/CRCvsSTAD_NC/greater/CRC_vs_STAD.csv",header = TRUE,row.names = 1)
  
  CRC_highexpression_vs_NC <- CRCvsNC[CRCvsNC$p.value<0.05&CRCvsNC$log2FC>0.9,]
  CRC_highexpression_vs_STAD <- CRCvsSTAD[CRCvsSTAD$p.value<0.05&CRCvsSTAD$log2FC>0.9,]
  
  CRC <- intersect(CRC_highexpression_vs_NC$gene,CRC_highexpression_vs_STAD$gene)
  
  STADvsCRC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/STADvsCRC_NC/greater/STAD_vs_CRC.csv",header = TRUE,row.names = 1)
  STADvsNC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/STADvsCRC_NC/greater/STAD_vs_NC.csv",header = TRUE,row.names = 1)
  
  STAD_highexpression_vs_CRC <- STADvsCRC[STADvsCRC$p.value<0.05&STADvsCRC$log2FC>0.65,]
  STAD_highexpression_vs_NC <- STADvsNC[STADvsNC$p.value<0.05&STADvsNC$log2FC>0.65,]
  
  STAD <- intersect(STAD_highexpression_vs_CRC$gene,STAD_highexpression_vs_NC$gene)
  
  NCvsCRC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/NCvsCRC_STAD/greater/NC_vs_CRC.csv",header = TRUE,row.names = 1)
  NCvsSTAD <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/NCvsCRC_STAD/greater/NC_vs_STAD.csv",header = TRUE,row.names = 1)
  
  
  NC_highexpression_vs_CRC <- NCvsCRC[NCvsCRC$p.value<0.05&NCvsCRC$log2FC>1.4,]
  NC_highexpression_vs_STAD <- NCvsSTAD[NCvsSTAD$p.value<0.05&NCvsSTAD$log2FC>1.4,]
  
  NC <- intersect(NC_highexpression_vs_CRC$gene,NC_highexpression_vs_STAD$gene)
  
  CRC <- as.data.frame(CRC)
  colnames(CRC) <- "gene"
  STAD <- as.data.frame(STAD)
  colnames(STAD) <- "gene"
  NC <- as.data.frame(NC)
  colnames(NC) <- "gene"
  
  
  highexpression_genes <- rbind(CRC,STAD,NC)
}
#less
{
  CRCvsNC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/CRCvsSTAD_NC/less/CRC_vs_NC.csv",header = TRUE,row.names = 1)
  CRCvsSTAD <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/CRCvsSTAD_NC/less/CRC_vs_STAD.csv",header = TRUE,row.names = 1)
  
  CRC_highexpression_vs_NC <- CRCvsNC[CRCvsNC$p.value<0.05&CRCvsNC$log2FC< -0.95,]
  CRC_highexpression_vs_STAD <- CRCvsSTAD[CRCvsSTAD$p.value<0.05&CRCvsSTAD$log2FC< -0.95,]
  
  CRC <- intersect(CRC_highexpression_vs_NC$gene,CRC_highexpression_vs_STAD$gene)
  
  STADvsCRC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/STADvsCRC_NC/less/STAD_vs_CRC.csv",header = TRUE,row.names = 1)
  STADvsNC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/STADvsCRC_NC/less/STAD_vs_NC.csv",header = TRUE,row.names = 1)
  
  STAD_highexpression_vs_CRC <- STADvsCRC[STADvsCRC$p.value<0.05&STADvsCRC$log2FC< -0.7,]
  STAD_highexpression_vs_NC <- STADvsNC[STADvsNC$p.value<0.05&STADvsNC$log2FC< -0.7,]
  
  STAD <- intersect(STAD_highexpression_vs_CRC$gene,STAD_highexpression_vs_NC$gene)
  
  NCvsCRC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/NCvsCRC_STAD/less/NC_vs_CRC.csv",header = TRUE,row.names = 1)
  NCvsSTAD <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/NCvsCRC_STAD/less/NC_vs_STAD.csv",header = TRUE,row.names = 1)
  
  
  NC_highexpression_vs_CRC <- NCvsCRC[NCvsCRC$p.value<0.05&NCvsCRC$log2FC< -1.2,]
  NC_highexpression_vs_STAD <- NCvsSTAD[NCvsSTAD$p.value<0.05&NCvsSTAD$log2FC< -1.2,]
  
  NC <- intersect(NC_highexpression_vs_CRC$gene,NC_highexpression_vs_STAD$gene)
  
  CRC <- as.data.frame(CRC)
  colnames(CRC) <- "gene"
  STAD <- as.data.frame(STAD)
  colnames(STAD) <- "gene"
  NC <- as.data.frame(NC)
  colnames(NC) <- "gene"
  
  
  highexpression_genes <- rbind(CRC,STAD,NC)
}
#edgeR_glmlrt CRCvsSTADvsNC
{
CRCvsSTAD_NC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/edgeR/edgeR_glmlrt-CRCvsSTAD_NC.txt",header = TRUE,sep = "\t",row.names = 1)

CRC_highexpression <- CRCvsSTAD_NC[CRCvsSTAD_NC$padj<0.0001&CRCvsSTAD_NC$log2FoldChange>0,]

CRC_highexpression <- CRC_highexpression[grep("mRNA|lncRNA",rownames(CRC_highexpression)),]

CRC_highexpression <- head(CRC_highexpression[order(CRC_highexpression$log2FoldChange,decreasing = T),],10)

nrow(CRC_highexpression)

STADvsCRC_NC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/edgeR/edgeR_glmlrt-STADvsCRC_NC.txt",header = TRUE,sep = "\t",row.names = 1)

STAD_highexpression <- STADvsCRC_NC[STADvsCRC_NC$padj<0.0001&STADvsCRC_NC$log2FoldChange>0,]

STAD_highexpression <- STAD_highexpression[grep("mRNA|lncRNA",rownames(STAD_highexpression)),]

STAD_highexpression <- head(STAD_highexpression[order(STAD_highexpression$log2FoldChange,decreasing = T),],10)

nrow(STAD_highexpression)

NCvsCRC_STAD <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/edgeR/edgeR_glmlrt-NCvsCRC_STAD.txt",header = TRUE,sep = "\t",row.names = 1)

NC_highexpression <- NCvsCRC_STAD[NCvsCRC_STAD$padj<0.0001&NCvsCRC_STAD$log2FoldChange>0,]

NC_highexpression <- NC_highexpression[grep("mRNA|lncRNA",rownames(NC_highexpression)),]

NC_highexpression <- head(NC_highexpression[order(NC_highexpression$log2FoldChange,decreasing = T),],10)

nrow(NC_highexpression)

intersect(rownames(CRC_highexpression),rownames(STAD_highexpression))

CRC <- as.data.frame(rownames(CRC_highexpression))
colnames(CRC) <- "gene"
STAD <- as.data.frame(rownames(STAD_highexpression))
colnames(STAD) <- "gene"
NC <- as.data.frame(rownames(NC_highexpression))
colnames(NC) <- "gene"

highexpression_genes <- rbind(CRC,STAD,NC)
}


  
#heatmap
#CRC vs STAD vs NC
library(pheatmap)
{
  target <- highexpression_genes$gene
  type <- c("CRC","STAD","NC")
  counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/multiomics_TPM.txt",header = T, row.names = 1,sep="\t")
  
  df.type=as.data.frame(array(dim=c(1,1))) #设置数据框行列
  i=1
  while(i<=length(type)){
    df <- counts[,grep(type[i],colnames(counts),fixed=TRUE)]
    df.type <- cbind(df.type,df)
    i=i+1
  }
  df.type <- df.type[,-1]
  
  df.type <- df.type[,-grep("CRC.PKU.34|STAD.PKU.35",colnames(df.type))]
  
  # progress
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = length(target), clear = FALSE, width= 60) 
  df.target={}
  j=1
  while(j<=length(target)){
    df <- df.type[grep(target[j],rownames(df.type),fixed=TRUE),]
    df.target <- rbind(df.target,df)
    j=j+1
    pb$tick()
    Sys.sleep(1 / 100)
  }
  
  heatmap <- df.target
  
  #for inhouse processed datasets, necessary to extract ENSEMBL ID, for downloaded processed dataset, not necessary
  rownames(heatmap) <- as.character(lapply(strsplit(rownames(df.target),"\\|"),function(x) x[3])) #1 for ENSG at each line begin, 3 for pico(ENSG at the third col)
  # key type transvert
  #symbol <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(heatmap), columns=c("SYMBOL"), keytype="ENSEMBL")
  #symbol$ENSEMBL #check and find 77st row ensembl ID repeatly show up(for gut GSE133864) #rownames(heatmap) <- symbol[-78,]$SYMBOL
  #rownames(heatmap) <- symbol[-33,]$SYMBOL #check and find 32st row ensembl ID repeatly show up(for GSE89843)
  #rownames(heatmap) <- symbol[-76,]$SYMBOL #check and find 76st row ensembl ID repeatly show up(for GSE68086)
  #rownames(heatmap) <- symbol[-119,]$SYMBOL #check and find 119st row ensembl ID repeatly show up(for pico)
  #rownames(heatmap) <- symbol[-136,]$SYMBOL #check and find 119st row ensembl ID repeatly show up(for THCA)
  
  column_annotation = as.data.frame(array(dim=c(length(colnames(heatmap)),1))) 
  rownames(column_annotation) <- colnames(heatmap)
  #GSE68086 must manully revise
  #write.csv(colnames(heatmap),"heatmap_annotation.csv")
  #column_annotation <- read.csv("heatmap_annotation.csv",row.names = 1,header = F)
  column_annotation$V1 <- as.character(lapply(strsplit(colnames(heatmap),".",fixed = TRUE),function(x) x[1])) ##exoRBase "_",x[1]; pico ".",x[1];GSE89843, ".", x[2]
  colnames(column_annotation) <- c("Type")
  
  
  bk = unique(c(seq(-1,1, length=100)))
  pheatmap(
    heatmap,
    annotation_col = column_annotation,
    scale = "row",
    cluster_cols = FALSE,cluster_rows = FALSE,
    show_colnames= FALSE, 
    #cluster_cols = FALSE,
    breaks=bk,
    gaps_row = c(nrow(CRC),nrow(CRC)+nrow(STAD)),
    gaps_col = c(38,73),
    colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 10)
}

#CRC_STADvsNC
{
{
CRC_STADvsNC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/edgeR/edgeR_glmlrt-CRC_STADvsNC.txt",header = TRUE,sep = "\t",row.names = 1)

CRC_STAD_highexpression <- CRC_STADvsNC[CRC_STADvsNC$padj<0.05&CRC_STADvsNC$log2FoldChange>0,]

CRC_STAD_highexpression <- CRC_STAD_highexpression[grep("mRNA|lncRNA",rownames(CRC_STAD_highexpression)),]

CRC_STAD_highexpression <- head(CRC_STAD_highexpression[order(CRC_STAD_highexpression$log2FoldChange,decreasing = T),],50)

CRC_STAD_lowexpression <- CRC_STADvsNC[CRC_STADvsNC$padj<0.05&CRC_STADvsNC$log2FoldChange<0,]

CRC_STAD_lowexpression <- CRC_STAD_lowexpression[grep("mRNA|lncRNA",rownames(CRC_STAD_lowexpression)),]

CRC_STAD_lowexpression <- head(CRC_STAD_lowexpression[order(CRC_STAD_lowexpression$log2FoldChange,decreasing = T),],50)

nrow(CRC_STAD_highexpression)
nrow(CRC_STAD_lowexpression)

CRC_STAD <- as.data.frame(rownames(CRC_STAD_highexpression))
colnames(CRC_STAD) <- "gene"
NC <- as.data.frame(rownames(CRC_STAD_lowexpression))
colnames(NC) <- "gene"

highexpression_genes <- rbind(CRC_STAD,NC)
}

{
  target <- highexpression_genes$gene
  type <- c("CRC","STAD","NC")
  counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/multiomics_TPM.txt",header = T, row.names = 1,sep="\t")
  
  df.type=as.data.frame(array(dim=c(1,1))) #设置数据框行列
  i=1
  while(i<=length(type)){
    df <- counts[,grep(type[i],colnames(counts),fixed=TRUE)]
    df.type <- cbind(df.type,df)
    i=i+1
  }
  df.type <- df.type[,-1]
  
  df.type <- df.type[,-grep("CRC.PKU.34|STAD.PKU.35",colnames(df.type))]
  
  # progress
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = length(target), clear = FALSE, width= 60) 
  df.target={}
  j=1
  while(j<=length(target)){
    df <- df.type[grep(target[j],rownames(df.type),fixed=TRUE),]
    df.target <- rbind(df.target,df)
    j=j+1
    pb$tick()
    Sys.sleep(1 / 100)
  }
  
  heatmap <- df.target
  
  #for inhouse processed datasets, necessary to extract ENSEMBL ID, for downloaded processed dataset, not necessary
  rownames(heatmap) <- as.character(lapply(strsplit(rownames(df.target),"\\|"),function(x) x[3])) #1 for ENSG at each line begin, 3 for pico(ENSG at the third col)
  # key type transvert
  #symbol <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(heatmap), columns=c("SYMBOL"), keytype="ENSEMBL")
  #symbol$ENSEMBL #check and find 77st row ensembl ID repeatly show up(for gut GSE133864) #rownames(heatmap) <- symbol[-78,]$SYMBOL
  #rownames(heatmap) <- symbol[-33,]$SYMBOL #check and find 32st row ensembl ID repeatly show up(for GSE89843)
  #rownames(heatmap) <- symbol[-76,]$SYMBOL #check and find 76st row ensembl ID repeatly show up(for GSE68086)
  #rownames(heatmap) <- symbol[-119,]$SYMBOL #check and find 119st row ensembl ID repeatly show up(for pico)
  #rownames(heatmap) <- symbol[-136,]$SYMBOL #check and find 119st row ensembl ID repeatly show up(for THCA)
  
  column_annotation = as.data.frame(array(dim=c(length(colnames(heatmap)),1))) 
  rownames(column_annotation) <- colnames(heatmap)
  #GSE68086 must manully revise
  #write.csv(colnames(heatmap),"heatmap_annotation.csv")
  #column_annotation <- read.csv("heatmap_annotation.csv",row.names = 1,header = F)
  column_annotation$V1 <- as.character(lapply(strsplit(colnames(heatmap),".",fixed = TRUE),function(x) x[1])) ##exoRBase "_",x[1]; pico ".",x[1];GSE89843, ".", x[2]
  colnames(column_annotation) <- c("Type")
  
  
  bk = unique(c(seq(-1,1, length=100)))
  pheatmap(
    heatmap,
    annotation_col = column_annotation,
    scale = "row",
    cluster_cols = FALSE,cluster_rows = FALSE,
    show_colnames= FALSE, 
    #cluster_cols = FALSE,
    breaks=bk,
    gaps_row = c(nrow(CRC_STAD)),
    gaps_col = c(73),
    colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 6)
}
}

#CRCvsSTAD
{
  {
    CRCvsSTAD <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/edgeR/edgeR_glmlrt-CRCvsSTAD.txt",header = TRUE,sep = "\t",row.names = 1)
    
    CRCvsSTAD_highexpression <- CRCvsSTAD[CRCvsSTAD$padj<0.05&CRCvsSTAD$log2FoldChange>0,]
    
    CRCvsSTAD_highexpression <- CRCvsSTAD_highexpression[grep("mRNA|lncRNA",rownames(CRCvsSTAD_highexpression)),]
    
    CRCvsSTAD_highexpression <- head(CRCvsSTAD_highexpression[order(CRCvsSTAD_highexpression$log2FoldChange,decreasing = T),],50)
    
    CRCvsSTAD_lowexpression <- CRCvsSTAD[CRCvsSTAD$padj<0.05&CRCvsSTAD$log2FoldChange<0,]
    
    CRCvsSTAD_lowexpression <- CRCvsSTAD_lowexpression[grep("mRNA|lncRNA",rownames(CRCvsSTAD_lowexpression)),]
    
    CRCvsSTAD_lowexpression <- head(CRCvsSTAD_lowexpression[order(CRCvsSTAD_lowexpression$log2FoldChange,decreasing = T),],50)
    
    nrow(CRCvsSTAD_highexpression)
    nrow(CRCvsSTAD_lowexpression)
    
    CRC <- as.data.frame(rownames(CRCvsSTAD_highexpression))
    colnames(CRC) <- "gene"
    STAD <- as.data.frame(rownames(CRCvsSTAD_lowexpression))
    colnames(STAD) <- "gene"
    
    highexpression_genes <- rbind(CRC,STAD)
  }
  
  {
    target <- highexpression_genes$gene
    type <- c("CRC","STAD")
    counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/multiomics_TPM.txt",header = T, row.names = 1,sep="\t")
    
    df.type=as.data.frame(array(dim=c(1,1))) #设置数据框行列
    i=1
    while(i<=length(type)){
      df <- counts[,grep(type[i],colnames(counts),fixed=TRUE)]
      df.type <- cbind(df.type,df)
      i=i+1
    }
    df.type <- df.type[,-1]
    
    df.type <- df.type[,-grep("CRC.PKU.34|STAD.PKU.35",colnames(df.type))]
    
    # progress
    pb <- progress_bar$new(
      format = "  Processing [:bar] :percent eta: :eta",
      total = length(target), clear = FALSE, width= 60) 
    df.target={}
    j=1
    while(j<=length(target)){
      df <- df.type[grep(target[j],rownames(df.type),fixed=TRUE),]
      df.target <- rbind(df.target,df)
      j=j+1
      pb$tick()
      Sys.sleep(1 / 100)
    }
    
    heatmap <- df.target
    
    #for inhouse processed datasets, necessary to extract ENSEMBL ID, for downloaded processed dataset, not necessary
    rownames(heatmap) <- as.character(lapply(strsplit(rownames(df.target),"\\|"),function(x) x[3])) #1 for ENSG at each line begin, 3 for pico(ENSG at the third col)
    # key type transvert
    #symbol <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(heatmap), columns=c("SYMBOL"), keytype="ENSEMBL")
    #symbol$ENSEMBL #check and find 77st row ensembl ID repeatly show up(for gut GSE133864) #rownames(heatmap) <- symbol[-78,]$SYMBOL
    #rownames(heatmap) <- symbol[-33,]$SYMBOL #check and find 32st row ensembl ID repeatly show up(for GSE89843)
    #rownames(heatmap) <- symbol[-76,]$SYMBOL #check and find 76st row ensembl ID repeatly show up(for GSE68086)
    #rownames(heatmap) <- symbol[-119,]$SYMBOL #check and find 119st row ensembl ID repeatly show up(for pico)
    #rownames(heatmap) <- symbol[-136,]$SYMBOL #check and find 119st row ensembl ID repeatly show up(for THCA)
    
    column_annotation = as.data.frame(array(dim=c(length(colnames(heatmap)),1))) 
    rownames(column_annotation) <- colnames(heatmap)
    #GSE68086 must manully revise
    #write.csv(colnames(heatmap),"heatmap_annotation.csv")
    #column_annotation <- read.csv("heatmap_annotation.csv",row.names = 1,header = F)
    column_annotation$V1 <- as.character(lapply(strsplit(colnames(heatmap),".",fixed = TRUE),function(x) x[1])) ##exoRBase "_",x[1]; pico ".",x[1];GSE89843, ".", x[2]
    colnames(column_annotation) <- c("Type")
    
    
    bk = unique(c(seq(-1,1, length=100)))
    pheatmap(
      heatmap,
      annotation_col = column_annotation,
      scale = "row",
      cluster_cols = FALSE,cluster_rows = FALSE,
      show_colnames= FALSE, 
      #cluster_cols = FALSE,
      breaks=bk,
      gaps_row = c(nrow(CRC)),
      gaps_col = c(38),
      colorRampPalette(c("blue","white","red"))(100),
      fontsize_row = 6)
  }
}

##KEGG
{
STAD_NC_KEGG <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/KEGGandGSEA/CRC_STADvsNC_KEGG.txt",header=T,sep="\t",quote = "")

STAD_NC_KEGG$Term_short <- as.character(lapply(strsplit(as.character(STAD_NC_KEGG$Term),":",fixed = TRUE),function(x) x[2]))

KEGG <- head(STAD_NC_KEGG[order(STAD_NC_KEGG$PValue,decreasing = FALSE),],10)

KEGG$Term_short <- factor(KEGG$Term_short,levels=KEGG[order(KEGG$PValue,decreasing = TRUE),]$Term_short)

ggplot(KEGG,aes(x=Fold.Enrichment,y=Term_short))+
  geom_point(aes(size=Count,color=-1*log10(PValue)))+
  scale_colour_gradient(low="green",high="red")+
  labs(color=expression(-log[10](P.value)),
       size="Gene number",
       x="Fold enrichment")+
  theme_bw()+
  theme(axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
        axis.title.x = element_text(size=rel(1.3)),
        axis.title.y = element_blank())
}

##insert length
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/level_3_features/Insert_length/miso/")
  library(reshape2)
  library(ggplot2)
  library(extrafont)
  library(plyr)
  ###insert length by miso
  
  plot <- read.csv("miso_Insert_length_summary_pico.txt",sep = "\t",header = T)
  plot <- plot[-1,]
  
  ##length density
  {
    colnames(plot)
    #c(rep("Gastric Cancer",15),rep("Healthy Donor",15))
    #colnames(plot) = c(rep("STAD",15),rep("Healthy_Donor",15))
    #colnames(plot) = c(rep("sample",30))
    plot.m <- melt(plot,na.rm = TRUE)
    head(plot.m)
    colnames(plot.m) = c("Sample","Insert Length")
    p<-ggplot(plot.m, aes(x = plot.m$`Insert Length`))+
      theme_bw(base_size = 12, base_family = "Arial")+
      theme(#legend.position="right",
        legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        axis.text.x = element_text(face="bold", color="black", size=20),
        axis.text.y = element_text(face="bold",  color="black", size=20),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      labs(x="Fragment Size(nt)",y="Density",title="", face="bold")
    
    p+geom_density(aes(fill = Sample),alpha = 0.4,size = 0)+stat_density(adjust = 1)
    p+geom_histogram(bin=10)
  }
  
  {
    # freqency line plot
    summary(plot[,1]) #看数据的分布区间
    # set bin, min, max
    bin=1
    min=1
    max=1000
    m <- seq(min-1,max,by=bin)#设置一个区间范围
    barplot(table(cut(plot[,1],m)))#简单作图查看分布
    
    
    # 计算frequency
    freq = {}
    col <- colnames(plot) 
    i=1
    while(i<=length(col)){
      #prop.table(table(cut(plot$col[i],m)))*100#计算各个区间的频率
      freq.tmp <- as.data.frame(prop.table(table(cut(plot[,col[i]],m)))*100)
      freq <- cbind(freq,freq.tmp$Freq)
      i=i+1
    }
    # make frequency dataframe
    freq <- as.data.frame(freq)
    colnames(freq) <- colnames(plot)
    rownames(freq) <- as.data.frame(prop.table(table(cut(plot[,1],m))))[,1]
    
    # summary group frequency
    summary= {}
    summary$mean <- apply(freq, 1, mean,na.rm=T)
    #summary$mean_STAD <- apply(freq[,grep("STAD",colnames(freq))], 1, mean,na.rm=T)
    #summary$mean_HCC <- apply(freq[,grep("HCC",colnames(freq))], 1, mean,na.rm=T)
    #summary$mean_CRC <- apply(freq[,grep("CRC",colnames(freq))], 1, mean,na.rm=T)
    #summary$mean_ESCA <- apply(freq[,grep("ESCA",colnames(freq))], 1, mean,na.rm=T)
    #summary$mean_LUAD <- apply(freq[,grep("LUAD",colnames(freq))], 1, mean,na.rm=T)
    #summary$mean_NC <- apply(freq[,grep("NC",colnames(freq))], 1, mean,na.rm=T)
    summary.df <- as.data.frame(summary)
    summary.df$length <- seq(from=min-1,to=max-bin,by=bin)
    write.csv(summary.df,"miso_insert_pico.csv")
    summary.df <- read.csv("miso_insert_compare.csv",header=T)
    
    #plot
    summary.df$Dataset <- factor(summary.df$Dataset,levels=c("DNA","pico","GSE133684","exoRBase"))
    p <- ggplot(summary.df, aes(x=length,y=mean,color = Dataset)) + geom_line(size = 2) + theme_bw(base_size = 12, base_family = "Arial")+
      geom_vline(aes(xintercept = 167),linetype="dashed")+
      #scale_color_brewer(palette="Blues") + 
      scale_color_manual(values=c("#FFCC11","#5CACEE","#104E8B","#687C97"))+
      scale_x_continuous(limits=c(0, 1000),breaks = c(0,100,200,400,600,800,1000))+
      #xlim(0,500)+
      #scale_y_continuous(limits=c(0, 10),breaks = c(0.0,2.5,5,7.5,10))+
      theme(#legend.position="right",
        legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=30),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=30),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_text(face="bold", color="black", size=30),
        axis.text.y = element_text(face="bold",  color="black", size=30),
        axis.title.x = element_text(face="bold", color="black", size=36),
        axis.title.y = element_text(face="bold",color="black", size=36))+
      labs(x="Insert Size (nt)",y="Freuqency (%)",title="", face="bold")
    
    ggsave("miso_insert_length_compare.pdf",p)
  }
}




#RP-RNA multi features (and cancer genes/Neutral genes )
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/RP-RNA features/feature matrix/")
  #get targeted genes
  {
  RP_RNA <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/RP-RNA features/genes ploted/genes_RP-mRNA_COSMIC_Neutral.csv",header = T)
  
  #Alt.promoter
  counts <- read.csv("Alt.promoter_TPM-by-promoter_multiomics_paired.txt",sep="\t",header = T, row.names = 1)
  #APA
  counts <- read.csv("APA_PDUI_CRC)STADvsNC.txt",sep="\t",header = T, row.names = 1)
  #chimeric
  counts <- read.csv("chimeric_counts.txt",sep="\t",header = T, row.names = 1)
  ##splicing
  #A3SS
  counts <- read.csv("splicing/A3SS_JC_inc_level.txt",sep="\t",header = T, row.names = 1)
  #A5SS
  counts <- read.csv("splicing/A5SS_JC_inc_level.txt",sep="\t",header = T, row.names = 1)
  #MXE
  counts <- read.csv("splicing/MXE_JC_inc_level.txt",sep="\t",header = T, row.names = 1)
  #RI
  counts <- read.csv("splicing/RI_JC_inc_level.txt",sep="\t",header = T, row.names = 1)
  #SE
  counts <- read.csv("splicing/SE_JC_inc_level.txt",sep="\t",header = T, row.names = 1)
  
  #get 78 RP RNAs and otehr cancer or neutral genes
  j=1
  pathway_gene_count={}
  while(j<=nrow(RP_RNA)){
    target <- RP_RNA[j,2]
    gene_symbol <- RP_RNA[j,1]
    if(length(grep(target,rownames(counts)))==0) {
      #print(paste0("No ",target," in this dataset."))
      temp <- as.data.frame(array(,dim=c(1,ncol(counts))))
      temp[1,] <- 0
      rownames(temp) <- gene_symbol
      colnames(temp) <- colnames(counts)
      pathway_gene_count <- rbind(pathway_gene_count,temp)
      j=j+1
    } else {
      #temp <- counts[which(rownames(counts)==target),]#for gene symbol
      temp <- counts[grep(target,rownames(counts),fixed=TRUE),]  #for ensg
      #temp <- counts[grep(target,rownames(counts),fixed=TRUE),] 
      #rownames(temp) <- gene_symbol
      pathway_gene_count <- rbind(pathway_gene_count,temp)
      j=j+1
    }
  }
  write.csv(pathway_gene_count,"../RP-mRNA_COSMIC_Neutral/SE_RP-mRNA_COSMIC_Neutral.csv")
  }
  
  {
  splicing_A3SS <- read.csv("../RP-mRNA_COSMIC_Neutral/A3SS_RP-mRNA_COSMIC_Neutral.csv",header = T)
  splicing_A3SS$type <- "splicing_A3SS"
  splicing_A5SS <- read.csv("../RP-mRNA_COSMIC_Neutral/A5SS_RP-mRNA_COSMIC_Neutral.csv",header = T)
  splicing_A5SS$type <- "splicing_A5SS"
  splicing_MXE <- read.csv("../RP-mRNA_COSMIC_Neutral/MXE_RP-mRNA_COSMIC_Neutral.csv",header = T)
  splicing_MXE$type <- "splicing_MXE"
  splicing_RI <- read.csv("../RP-mRNA_COSMIC_Neutral/RI_RP-mRNA_COSMIC_Neutral.csv",header = T)
  splicing_RI$type <- "splicing_RI"
  splicing_SE <- read.csv("../RP-mRNA_COSMIC_Neutral/SE_RP-mRNA_COSMIC_Neutral.csv",header = T)
  splicing_SE$type <- "splicing_SE"
  
  alt.promoter <- read.csv("../RP-mRNA_COSMIC_Neutral/Alt.promoter_RP-mRNA_COSMIC_Neutral.csv",header = T)
  alt.promoter$type <- "alt.promoter"
  APA <- read.csv("../RP-mRNA_COSMIC_Neutral/APA_RP-mRNA_COSMIC_Neutral.csv",header = T)
  APA$type <- "APA"
  #SNP <- read.csv("./SNP/SNP_genelevel_count.csv",header = T)
  #SNP$type <- "SNP"
  }
  large_matrix <- plyr::rbind.fill(splicing_A3SS,splicing_A5SS,splicing_MXE,splicing_RI,splicing_SE,alt.promoter,APA)
  
  large_matrix$rownames <- paste(large_matrix$X,"-",large_matrix$type)
  
  row_annotation = as.data.frame(array(dim=c(length(large_matrix$rownames),1))) 
  rownames(row_annotation) <- large_matrix$rownames
  row_annotation$V1 <- large_matrix$type
  colnames(row_annotation) <- c("Type")
  
  
  
  
  #overall
  plot <- large_matrix[,2:98]
  rownames(plot) <- large_matrix$rownames
  
  #single feature
  plot <- large_matrix[grep("APA",large_matrix$type),2:98]
  rownames(plot) <- large_matrix[grep("APA",large_matrix$type),1]
  
  
 
  plot[is.na(plot)] <- 0
  bk = unique(c(seq(-0.5,0.5, length=100)))
  
  col_annotation <- as.data.frame(array(dim=c(length(colnames(plot)),1)))
  rownames(col_annotation) <- colnames(plot)
  col_annotation$V1 <-  as.character(lapply(strsplit(colnames(plot),"\\."),function(x) x[1]))
  colnames(col_annotation) <- c("Type")
  col_annotation$Type <- factor(col_annotation$Type, levels = c("CRC","NC","STAD"))
  
  pheatmap(plot,
           breaks=bk,
           colorRampPalette(c("blue","white","red"))(100),
           annotation_row = row_annotation,
           #annotation_col = col_annotation,
           cluster_cols = F,cluster_rows = F,
           show_rownames=FALSE,
           show_colnames=FALSE,
           scale = "row",
           )
  
  
 
  
  
}

##waterfall plot by genVisR
{
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenVisR")

library("GenVisR")
set.seed(426)

# Plot with the MAF file type specified (default) The mainRecurCutoff parameter
# is described in the next section
waterfall(brcaMAF, fileType = "MAF", mainRecurCutoff = 0.05)
}

##waterfall plot by maftools
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("maftools")
  
  library(maftools)
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/SNP/annotation_filteredPBMC/")
  var.annovar.maf = annovarToMaf(annovar = "CRC_sample_CRC_COSMIC_cancergenes_filtered.txt", 
                                 Center = 'NA', 
                                 refBuild = 'hg38', 
                                 tsbCol = 'Tumor_Sample_Barcode', 
                                 table = 'refGene',
                                 sep = "\t")
  write.table(var.annovar.maf,file="var_annovar.maf",quote= F,sep="\t",row.names=F)
  var_maf = read.maf(maf ="var_annovar.maf")
  
  var_maf@clinical.data$Type <- as.character(lapply(strsplit(as.character(var_maf@clinical.data$Tumor_Sample_Barcode),"-"),function(x) x[1]))
  
  var_maf@clinical.data
  
  plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median')
  
  oncoplot(maf = var_maf, top = 25 ,showTumorSampleBarcodes = F,clinicalFeatures = "Type",colbar_pathway = FALSE,
           sortByAnnotation = TRUE,draw_titv = TRUE)
  
  laml.titv = titv(maf = var_maf, plot = FALSE, useSyn = TRUE)
  plotTiTv(res = laml.titv)
  
  somaticInteractions(maf = var_maf, top = 10, pvalue = c(0.05, 0.1))
  
  #lollipop plot for APC
  lollipopPlot(
    maf = var_maf,
    gene = 'APC',
    AACol = 'aaChange',
    showMutationRate = TRUE,
    #labelPos = "all",
    #refSeqID = "NM_000546",
    printCount = TRUE,
    showDomainLabel = TRUE
  )
  
  #oncogene pathway
  OncogenicPathways(maf = var_maf)
  PlotOncogenicPathways(maf = var_maf,pathways = "WNT")
  
  #
  rainfallPlot(maf = var_maf, detectChangePoints = TRUE, pointSize = 0.4)
}


#gene expression QC
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/QC")
TPM <- read.csv("TPM_intron_spanning_noMTRNA_passed_filter_noMIX.txt",header =TRUE, row.names =1,sep = "\t")

View(TPM)

TPM$rowsum <- rowSums(TPM)

sorted_TPM <- TPM[order(TPM$rowsum,decreasing = T),]

write.csv(sorted_TPM[1:15,],"TPM_top15.csv")




#Clustering
{
  #ConsensusClusterPlus
{
  library(ConsensusClusterPlus)
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/04.Clustering")
  title <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/04.Clustering/pam_byCC"
  DE_CRC_vs_NC <- read.csv("TPM_CRCvsNC_DE.txt",header = TRUE, row.names = 1,sep ="\t")
  results = ConsensusClusterPlus(as.matrix(DE_CRC_vs_NC),maxK=6,reps=50,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
  
  TPM_CRC <- read.csv("TPM_CRC_passedQC.txt",header = TRUE, row.names = 1,sep ="\t")
  mads=apply(TPM_CRC,1,mad)
  d=TPM_CRC[rev(order(mads))[1:5000],]
  d=TPM_CRC[rev(order(mads)),]
  d = sweep(d,1, apply(d,1,median,na.rm=T))
  results = ConsensusClusterPlus(as.matrix(d),maxK=6,reps=1000,pItem=0.8,pFeature=1,title=title,clusterAlg="pam",distance="pearson",seed=1262118388.71279,plot="png")
  results[[4]][["consensusMatrix"]][1:5,1:5]
  results[[4]][["consensusClass"]]
  
  
  col_anno <- as.data.frame(results[[4]][["consensusClass"]])
  colnames(col_anno) <- "Consensus Cluster"
  col_anno <- col_anno[order(col_anno$`Consensus Cluster`),,drop=FALSE]
  d <- d[,rownames(col_anno)]
  bk = unique(c(seq(-1,1, length=100)))
  pheatmap(d,
           breaks=bk,
           #scale = "row",
           show_colnames= FALSE,show_rownames= FALSE,
           cluster_rows = TRUE,cluster_cols = FALSE,
           annotation_col = col_anno)

  }

  #CMSclassifier
{
  install.packages("remotes")
  remotes::install_github("Sage-Bionetworks/CMSclassifier")
  library(CMSclassifier)
  sampleData <- read.table(d, sep="\t",header = TRUE,row.names = 1,check.names=FALSE)
  
  Rfcms <- CMSclassifier::classifyCMS(DE_CRC_vs_NC,method="RF")[[3]]
}

  #cola
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("cola")
  
  library(cola)
  
  mat = adjust_matrix(TPM_CRC)  # optional
  rl = run_all_consensus_partition_methods(mat, mc.cores = 1)
  cola_report(rl, output_dir = "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Clustering", mc.cores = 3)
}
}

##machine learning
{
library(caret)
library(randomForest)
library(pROC)
library(PRROC)
library(dplyr)
library(plyr)
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/")
outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20210722_multiomics/merged2/STADvsNC/"
positive <- "STAD"
negative <- "NC"
featureNumber <- 200
#trasnform data matrix, col are features and row are samples, then add class/label at last column

#single feature matrix input
{
#raw_matrix <- read.csv("Splicing_ML.txt",sep = "\t",header = TRUE, row.names = 1)
#full_matrix <- as.data.frame(t(raw_matrix))
}

#merged feature matrix input
{
  full_matrix <- merged_selected_feature_matrix #merged matrix
}

full_matrix$label <- unlist(lapply(strsplit(rownames(full_matrix),".",fixed = TRUE),function(x) x[1]))

positive_matrix <- full_matrix[which(full_matrix$label==positive),]
negative_matrix <- full_matrix[which(full_matrix$label==negative),]
matrix <- rbind(positive_matrix,negative_matrix)
#numeric matrix
matrix_nolabel <- matrix[,-which(colnames(matrix)=="label")]
matrix_nolabel[is.na(matrix_nolabel)] <- 0
matrix_nolabel_numeric <- mutate_all(matrix_nolabel, function(x) as.numeric(as.character(x)))
rownames(matrix_nolabel_numeric) <- rownames(matrix_nolabel)

labels <- as.factor(matrix$label)

#feature selection top200
rfFuncs$summary <- twoClassSummary

set.seed(123)
subsetSizes <- c(featureNumber)
seeds <- vector(mode = "list", length = nrow(matrix_nolabel_numeric)+1 )
for(i in 1:(nrow(matrix_nolabel_numeric))) seeds[[i]] <- sample.int(1000, length(subsetSizes) + 1)
seeds[[(nrow(matrix_nolabel_numeric)+1)]] <- sample.int(1000, 1)

set.seed(1)
rfectrl <- rfeControl(functions=rfFuncs,
                      seeds = seeds,
                      saveDetails = TRUE,
                      verbose = TRUE,
                      method="LOOCV")

rfe.results <- rfe(matrix_nolabel_numeric,labels,
                   sizes = c(featureNumber),
                   rfeControl = rfectrl,
                   metric = "ROC")

#get feature importance
y <- rfe.results$variables
finalImp <- ddply(y[, c("Overall", "var")], .(var), function(x) mean(x$Overall,na.rm = TRUE))
names(finalImp)[2] <- "Overall"
finalImp <- finalImp[order(finalImp$Overall, decreasing = TRUE),]

if(rfe.results$optsize==featureNumber) {
  print(paste0(featureNumber," is RFE suggested. Using top"))
  feature_selected <- predictors(rfe.results)
} else if(featureNumber<=ncol(matrix_nolabel_numeric)){
  print(paste0(featureNumber," is not suggested. Using top",featureNumber," as final feature."))
  feature_selected <- as.character(finalImp$var[1:featureNumber])
} else {
  total_feature_number <- ncol(matrix_nolabel_numeric)
  print(paste0(featureNumber," is larger than total feature number. Using total feature number:",total_feature_number," as final feature."))
  feature_selected <- as.character(finalImp$var[1:total_feature_number])
}


matrix_nolabel_selected <- matrix_nolabel_numeric[,feature_selected]
write.table(matrix_nolabel_selected,paste0(outdir,"matrix_nolabel_selected.txt"),quote = FALSE,sep="\t")
write.table(finalImp,paste0(outdir,"feature_importance.txt"),quote = FALSE,sep="\t")

i=1
prob_final <- {}
while(i<=nrow(matrix)) {
train_labels <- as.factor(matrix[-i,]$label)
#fix model parameters, such as mtry, ntree
tune_mtry <- tuneRF(matrix_nolabel_selected[-i,],train_labels) #tune mtry
tune_mtry <- as.data.frame(tune_mtry)
mtry_best <- tune_mtry[order(tune_mtry$OOBError,decreasing=FALSE),][1,1]

#train RF model use n_sample-1
RF_model <- randomForest(matrix_nolabel_selected[-i,],train_labels, mtry = mtry_best)

#predict leaved 1 sample by trained model 
test_labels <- factor(matrix[i,]$label,level=c(positive,negative))
predict_prob <- predict(RF_model,newdata = matrix_nolabel_selected[i,], type = "prob")
prob_temp <- as.data.frame(predict_prob)
prob_final <- rbind(prob_final,prob_temp)
i=i+1
}
write.table(prob_final,paste0(outdir,"prob_final_LOO.txt"),quote = FALSE,sep="\t")


#report ROC
predicted <- prob_final
roc.curve <- roc(labels,predicted[,which(colnames(predicted)==positive)])
ci.auc(roc.curve,conf.level = 0.95)
record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
write.table(record,paste0(outdir,"matrices_record.txt"),quote = FALSE,sep="\t")

#plot
predicted$labels <- unlist(lapply(strsplit(rownames(predicted),".",fixed = TRUE),function(x) x[1]))
predicted$labels <- gsub(positive,1,predicted$labels)
predicted$labels <- gsub(negative,0,predicted$labels)
predicted$labels <- as.numeric(predicted$labels)

pdf(paste0(outdir,"AUROC.pdf"),width=7.8,height=6)
roc_curve_prediction_df <- PRROC::roc.curve(scores.class0=predicted[,which(colnames(predicted)==positive)], weights.class0 = predicted$labels,
                                            curve=TRUE,rand.compute = TRUE) # rand.compute 是0.5的随机分类器，为了画对角线的

par(cex.axis=2, cex.lab =1.6,font.axis =1.6, mar = c(5,5,2,3))
plot(roc_curve_prediction_df,rand.plot = TRUE,auc.main=FALSE,legend =4,lwd=5)
text(1.3, 0.2, paste('AUROC =', ' ',round(roc_curve_prediction_df$auc, digits = 4), sep =""), cex=2, font=2, col='black')
dev.off()

pdf(paste0(outdir,"AUPR.pdf"),width=7.8,height=6)
pr_curve_prediction_df <- PRROC::pr.curve(scores.class0=predicted[,which(colnames(predicted)==positive)], weights.class0 = predicted$labels,
                                            curve=TRUE) # rand.compute 是0.5的随机分类器，为了画对角线的

par(cex.axis=2, cex.lab =1.6,font.axis =1.6, mar = c(5,5,2,3))
plot(pr_curve_prediction_df,rand.plot = TRUE,auc.main=FALSE,legend =4,lwd=5)
text(1.3, 0.2, paste('AUPR =', ' ',round(pr_curve_prediction_df$auc.integral, digits = 4), sep =""), cex=2, font=2, col='black')
dev.off()

pred <- predicted[,which(colnames(predicted)==positive)]
ref <- unlist(lapply(strsplit(rownames(predicted),".",fixed = TRUE),function(x) x[1]))

i=1
while(i<=length(pred)){
  if(pred[i]>=0.5){
    pred[i] <- positive
  } else {
    pred[i] <- negative
  }
  i=i+1
}

pred <- factor(pred,levels = c(positive,negative))
ref <- factor(ref,levels = c(positive,negative))

confusion_matrix <- confusionMatrix(pred,ref)
write.table(confusion_matrix,paste0(outdir,"confusion_matrix.txt"),quote = FALSE,sep="\t")

}

#machine learning merged model
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20210722_multiomics/")
  all_RNA_features <- c("Alternative promoter outliers"="Alt.promoter",
                    "Expression outliers"="Expression",
                    "Splice outliers"="Splicing_FDR0.05",
                    "Alternative polyadenyltion outliers"="APA",
                    "Chimeric RNA"="chimeric",
                    "RNA editing outliers"="EDIT",
                    "ASE outliers"="ASE",
                    "RNA SNP"="SNP")
  
  i=1
  Group <- "STADvsNC"
  merged_selected_feature_matrix <- {}
  merged_selected_feature_matrix <- as.data.frame(matrix(numeric(0),nrow=54)) #47 for CRCvsNC; 53 for CRCvsSTAD; 54 for STADvsNC
  
  while (i <= length(all_RNA_features)){
    selected_feature_matrix <- read.csv(paste0(all_RNA_features[i],"/",Group,"/matrix_nolabel_selected.txt"),header = TRUE,row.names = 1,sep="\t")
    colnames(selected_feature_matrix) <- paste(as.character(all_RNA_features[i]), colnames(selected_feature_matrix), sep = "|")
    merged_selected_feature_matrix <- cbind(merged_selected_feature_matrix,selected_feature_matrix)
    i=i+1
  }
}

#plot different variants AUROC
{
  performance <- c("Alternative promoter"=0.9638,
                   "RNA Expression"=0.9321,
                   "RNA Splicing"=0.9022,
                   "Alternative polyadenyltion"=0.8714,
                   "Chimeric RNA"=0.8795,
                   "RNA editing"=0.9049,
                   "Allele specific expression"=0.8886,
                   "RNA SNP"=0.9692,
                   "Merged"=0.9783)
  
  RNA_color <- c("Alternative promoter"="#E9C2A6",
                "RNA Expression"="#A5435C",
                "RNA Splicing"="#C5E3BF",
                "Alternative polyadenyltion"="#003F87",
                "Chimeric RNA"="#FF3D0D",
                "RNA editing"="#324F17",
                "Allele specific expression"="#87CEFF",
                "RNA SNP"="#333333",
                "Merged"="white")
  performance_barplot <- as.data.frame(performance)
  performance_barplot$name <- rownames(performance_barplot)
  performance_barplot$name <- factor(performance_barplot$name,
                                     levels = performance_barplot[order(performance_barplot$performance,decreasing = TRUE),"name"])
  ggplot(performance_barplot,aes(x=name,y=performance-0.85,fill=name))+
    geom_bar(stat = "identity",colour = "black")+
    geom_text(aes(x=name,y=performance-0.85+0.005,label=performance),size = 4,angle = 0)+
    scale_y_continuous(breaks = c(0,0.05,0.10,0.15),labels = c("0.85","0.90","0.95","1"),expand = c(0,0),limits = c(0,0.15))+
    scale_fill_manual(values = RNA_color)+
    xlab("")+
    ylab("LOOCV AUROC")+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(15,5,10,5),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.ticks.x = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_text(face="bold", color="black", size=20),
      axis.title.y = element_text(face="bold",color="black", size=20))
    
}

#plot multi-ROC curve
library(pROC)
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20210722_multiomics/merged2/")

model_CRCvsNC <- read.csv("CRCvsNC/prob_final_LOO.txt",header = TRUE,row.names = 1,sep="\t")
model_CRCvsNC$response <- unlist(lapply(strsplit(rownames(model_CRCvsNC),".",fixed=TRUE),function(x) x[1]))
model_CRCvsNC$response <- gsub("CRC",1,model_CRCvsNC$response)
model_CRCvsNC$response <- gsub("NC",0,model_CRCvsNC$response)
model_STADvsNC <- read.csv("STADvsNC/prob_final_LOO.txt",header = TRUE,row.names = 1,sep="\t")
model_STADvsNC$response <- unlist(lapply(strsplit(rownames(model_STADvsNC),".",fixed=TRUE),function(x) x[1]))
model_STADvsNC$response <- gsub("STAD",1,model_STADvsNC$response)
model_STADvsNC$response <- gsub("NC",0,model_STADvsNC$response)
model_CRCvsSTAD <- read.csv("CRCvsSTAD/prob_final_LOO.txt",header = TRUE,row.names = 1,sep="\t")
model_CRCvsSTAD$response <- unlist(lapply(strsplit(rownames(model_CRCvsSTAD),".",fixed=TRUE),function(x) x[1]))
model_CRCvsSTAD$response <- gsub("CRC",1,model_CRCvsSTAD$response)
model_CRCvsSTAD$response <- gsub("STAD",0,model_CRCvsSTAD$response)

par(pty="s")
roc(model_CRCvsNC$response,model_CRCvsNC$CRC,plot=TRUE,legacy.axes=TRUE,
    xlab="False Positive Percentage",ylab="True Positive Percentage",col="#377eb8", lwd=4,print.auc=TRUE,print.auc.y= 0.25,print.auc.x= 0.35)
plot.roc(model_CRCvsSTAD$response,model_CRCvsSTAD$CRC,
         col="#EE9A00", lwd=4, print.auc=TRUE, print.auc.y=0.15, print.auc.x= 0.35, add=TRUE)
plot.roc(model_STADvsNC$response,model_STADvsNC$STAD,
         col="#4daf4a", lwd=4, print.auc=TRUE, print.auc.y=0.05, print.auc.x= 0.35, add=TRUE)
legend("bottomright",legend=c("CRCvsNC: 0.978","CRCvsSTAD: 0.978","STADvsNC: 0.983"),col=c("#377eb8","#EE9A00","#4daf4a"),lwd=4)


#process matrix
test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/featurecounts_genome_rmdup.txt",header = TRUE, row.names = 1, sep = "\t")
sample_ids <- colnames(test)

test$ID <- rownames(test)
test$Gene_symbol <- unlist(lapply(strsplit(rownames(test),"|",fixed = TRUE),function(x) x[3]))
test$Biotype <- unlist(lapply(strsplit(rownames(test),"|",fixed = TRUE),function(x) x[2]))
test$Ensembl_ID <- unlist(lapply(strsplit(rownames(test),".",fixed = TRUE),function(x) x[1]))
test$Transcript_length <- unlist(lapply(strsplit(rownames(test),"|",fixed = TRUE),function(x) x[7]))

output <- test[,c("ID","Ensembl_ID","Gene_symbol","Biotype","Transcript_length",sample_ids)]
colnames(output) <- gsub(".","-",colnames(output),fixed=TRUE)
write.table(output,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/Expression_demo.txt",quote = FALSE,sep = "\t",row.names = F)
