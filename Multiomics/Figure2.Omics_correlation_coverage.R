#correlation between plasma at sample level(reproducibility)
{
  TPM_Expression <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/PCAWG_biomarker_candidate/Expression/Expression_DNA_RNA_paired_with_MT.txt",sep="\t",header = T, row.names = 1)
  TPM_CNV <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/PCAWG_biomarker_candidate/DNA-CNV-CPM/CNV-CPM_DNA_RNA_paired.txt",sep="\t",header = T, row.names = 1)
  TPM_Methylation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/PCAWG_biomarker_candidate/MeDIP/MeDIP_DNA_RNA_paired.txt",sep="\t",header = T, row.names = 1)
  
  { 
    colnames(TPM_Expression) <- gsub(".","-",fixed=TRUE,colnames(TPM_Expression))
    colnames(TPM_Expression) <- gsub("-pico","-expression",fixed=TRUE,colnames(TPM_Expression))
    rownames(TPM_Expression) <- as.character(lapply(strsplit(rownames(TPM_Expression),"|",fixed = TRUE), function(x) x[1]))
    
    colnames(TPM_CNV) <- gsub(".","-",fixed=TRUE,colnames(TPM_CNV))
    colnames(TPM_CNV) <- gsub("-wgs","-cnv",fixed=TRUE,colnames(TPM_CNV))
    rownames(TPM_CNV) <- as.character(lapply(strsplit(rownames(TPM_CNV),"|",fixed = TRUE), function(x) x[1]))
    
    colnames(TPM_Methylation) <- gsub(".","-",fixed=TRUE,colnames(TPM_Methylation))
    colnames(TPM_Methylation) <- gsub("-me","-methylation",fixed=TRUE,colnames(TPM_Methylation))
    rownames(TPM_Methylation) <- as.character(lapply(strsplit(rownames(TPM_Methylation),"|",fixed = TRUE), function(x) x[1]))
    
    omics_paired <- cbind(TPM_CNV,TPM_Methylation)
    omics_paired <- inner_join(rownames_to_column(omics_paired),rownames_to_column(TPM_Expression),by=c("rowname"="rowname"))
    rownames(omics_paired) <- omics_paired$rowname
    omics_paired <- omics_paired[,-which(colnames(omics_paired)=="rowname")]
    omics_paired_ensembl <- omics_paired[grep("ENSG",rownames(omics_paired)),]
    
    # cor and pvalue bubble plot
    library(corrplot)
    library(RColorBrewer)
    
    M <-cor(omics_paired_ensembl)
    library("Hmisc")
    col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                   "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                   "#4393C3", "#2166AC", "#053061")))
    # Insignificant correlation are crossed
    res2 <- rcorr(as.matrix(omics_paired_ensembl),type="pearson")
    corrplot(corr = res2$r,col = col2(200),tl.col="black",type="lower", order="original",tl.pos = "ld",#tl.cex=0.7,tl.srt = 45,
             p.mat = res2$P, sig.level = 0.01,insig = "blank")  # insig = "blank"
    
    #write.csv(res2$r,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired.csv")
    #average_r <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average.csv",header = TRUE, row.names = 1)
    write.csv(res2$P,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired-20211007-p.csv")
    write.csv(res2$r,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired-20211007.csv")
    
    all_r <- res2$r
    all_p <- res2$P
    corrplot(corr = all_r[grep("cnv|methylation",rownames(all_r)),grep("cnv|methylation",colnames(all_r))],col = col2(200),tl.col="white",type="lower", order="original",tl.pos = "ld",#tl.cex=0.7,tl.srt = 45,
             p.mat = all_p[grep("cnv|methylation",rownames(all_p)),grep("cnv|methylation",colnames(all_p))], sig.level = 0.01,insig = "blank",addgrid.col = NA)  # insig = "blank"
    
    corrplot(corr = all_r[grep("cnv|expression",rownames(all_r)),grep("cnv|expression",colnames(all_r))],col = col2(200),tl.col="white",type="lower", order="original",tl.pos = "ld",#tl.cex=0.7,tl.srt = 45,
             p.mat = all_p[grep("cnv|expression",rownames(all_p)),grep("cnv|expression",colnames(all_p))], sig.level = 0.01,insig = "blank",addgrid.col = NA)  # insig = "blank"
    
    corrplot(corr = all_r[grep("methylation|expression",rownames(all_r)),grep("methylation|expression",colnames(all_r))],col = col2(200),tl.col="white",type="lower", order="original",tl.pos = "ld",#tl.cex=0.7,tl.srt = 45,
             p.mat = all_p[grep("methylation|expression",rownames(all_p)),grep("methylation|expression",colnames(all_p))], sig.level = 0.01,insig = "blank",addgrid.col = NA)  # insig = "blank"
    
    CNV_r <- all_r[grep("cnv",rownames(all_r)),grep("cnv",colnames(all_r))]
    average_CNV_r <- (sum(CNV_r)-ncol(CNV_r))/(ncol(CNV_r)*nrow(CNV_r)-ncol(CNV_r))
    
    Methy_r <- all_r[grep("methylation",rownames(all_r)),grep("methylation",colnames(all_r))]
    average_Methy_r <- (sum(Methy_r)-ncol(Methy_r))/(ncol(Methy_r)*nrow(Methy_r)-ncol(Methy_r))
    
    Expression_r <- all_r[grep("expression",rownames(all_r)),grep("expression",colnames(all_r))]
    average_Expression_r <- (sum(Expression_r)-ncol(Expression_r))/(ncol(Expression_r)*nrow(Expression_r)-ncol(Expression_r))
    
    Expression_CNV_r <- all_r[grep("expression",rownames(all_r)),grep("cnv",colnames(all_r))]
    average_Expression_CNV_r <- (sum(Expression_CNV_r))/(ncol(Expression_CNV_r)*nrow(Expression_CNV_r))
    
    Expression_Methy_r <- all_r[grep("expression",rownames(all_r)),grep("methylation",colnames(all_r))]
    average_Expression_Methy_r <- (sum(Expression_Methy_r))/(ncol(Expression_Methy_r)*nrow(Expression_Methy_r))
    
    CNV_Methy_r <- all_r[grep("cnv",rownames(all_r)),grep("methylation",colnames(all_r))]
    average_CNV_Methy_r <- (sum(CNV_Methy_r))/(ncol(CNV_Methy_r)*nrow(CNV_Methy_r))
    
    average_r <- as.data.frame(matrix(numeric(0),ncol=3,nrow=3))
    
    colnames(average_r) <- c("CNV","Methylation","Expression")
    rownames(average_r) <- c("CNV","Methylation","Expression")
    
    average_r[which(rownames(average_r)=="CNV"),which(colnames(average_r)=="CNV")] <- average_CNV_r
    average_r[which(rownames(average_r)=="Methylation"),which(colnames(average_r)=="Methylation")] <- average_Methy_r
    average_r[which(rownames(average_r)=="Expression"),which(colnames(average_r)=="Expression")] <- average_Expression_r
    average_r[which(rownames(average_r)=="Expression"),which(colnames(average_r)=="CNV")] <- average_Expression_CNV_r
    average_r[which(rownames(average_r)=="CNV"),which(colnames(average_r)=="Expression")] <- average_Expression_CNV_r
    average_r[which(rownames(average_r)=="Expression"),which(colnames(average_r)=="Methylation")] <- average_Expression_Methy_r
    average_r[which(rownames(average_r)=="Methylation"),which(colnames(average_r)=="Expression")] <- average_Expression_Methy_r
    average_r[which(rownames(average_r)=="CNV"),which(colnames(average_r)=="Methylation")] <- average_CNV_Methy_r
    average_r[which(rownames(average_r)=="Methylation"),which(colnames(average_r)=="CNV")] <- average_CNV_Methy_r
    
    
    write.csv(average_r,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average-20211007.csv")
    corrplot(corr = as.matrix(average_r),method = c("color"),col = col2(200),outline = TRUE,tl.col="black",type="lower", order="original",tl.pos = "l",tl.cex=1.5,cl.cex = 1, number.cex = as.matrix(average_r))
    
    
    {
      all_r <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired-20211007-p.csv",header = TRUE,row.names = 1)
      colnames(all_r) <- gsub(".","-",fixed = TRUE,colnames(all_r))
      Expression_r <- all_r[grep("expression",rownames(all_r)),grep("expression",colnames(all_r))]
      Methy_r <- all_r[grep("methylation",rownames(all_r)),grep("methylation",colnames(all_r))]
      CNV_r <- all_r[grep("cnv",rownames(all_r)),grep("cnv",colnames(all_r))]
      Expression_Methy_r <- all_r[grep("expression",rownames(all_r)),grep("methylation",colnames(all_r))]
      Expression_CNV_r <- all_r[grep("expression",rownames(all_r)),grep("cnv",colnames(all_r))]
      Methy_CNV_r <- all_r[grep("methylation",rownames(all_r)),grep("cnv",colnames(all_r))]
      
      ######Expression triangle
      Expression_CRC_r <- Expression_r[grep("CRC",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_r <- (sum(Expression_CRC_r))/(ncol(Expression_CRC_r)*nrow(Expression_CRC_r))
      
      Expression_STAD_r <- Expression_r[grep("STAD",rownames(Expression_r)),grep("STAD",colnames(Expression_r))]
      average_Expression_STAD_r <- (sum(Expression_STAD_r))/(ncol(Expression_STAD_r)*nrow(Expression_STAD_r))
      
      Expression_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("NC",colnames(Expression_r))]
      average_Expression_NC_r <- (sum(Expression_NC_r))/(ncol(Expression_NC_r)*nrow(Expression_NC_r))
      
      Expression_CRC_STAD_r <- Expression_r[grep("STAD",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_STAD_r <- (sum(Expression_CRC_STAD_r))/(ncol(Expression_CRC_STAD_r)*nrow(Expression_CRC_STAD_r))
      
      Expression_CRC_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_NC_r <- (sum(Expression_CRC_NC_r))/(ncol(Expression_CRC_NC_r)*nrow(Expression_CRC_NC_r))
      
      Expression_STAD_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("STAD",colnames(Expression_r))]
      average_Expression_STAD_NC_r <- (sum(Expression_STAD_NC_r))/(ncol(Expression_STAD_NC_r)*nrow(Expression_STAD_NC_r))
      
      ########Methylation triangle
      Methy_CRC_r <- Methy_r[grep("CRC",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_r <- (sum(Methy_CRC_r))/(ncol(Methy_CRC_r)*nrow(Methy_CRC_r))
      
      Methy_STAD_r <- Methy_r[grep("STAD",rownames(Methy_r)),grep("STAD",colnames(Methy_r))]
      average_Methy_STAD_r <- (sum(Methy_STAD_r))/(ncol(Methy_STAD_r)*nrow(Methy_STAD_r))
      
      Methy_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("NC",colnames(Methy_r))]
      average_Methy_NC_r <- (sum(Methy_NC_r))/(ncol(Methy_NC_r)*nrow(Methy_NC_r))
      
      Methy_CRC_STAD_r <- Methy_r[grep("STAD",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_STAD_r <- (sum(Methy_CRC_STAD_r))/(ncol(Methy_CRC_STAD_r)*nrow(Methy_CRC_STAD_r))
      
      Methy_CRC_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_NC_r <- (sum(Methy_CRC_NC_r))/(ncol(Methy_CRC_NC_r)*nrow(Methy_CRC_NC_r))
      
      Methy_STAD_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("STAD",colnames(Methy_r))]
      average_Methy_STAD_NC_r <- (sum(Methy_STAD_NC_r))/(ncol(Methy_STAD_NC_r)*nrow(Methy_STAD_NC_r))
      
      ########CNV triangle
      CNV_CRC_r <- CNV_r[grep("CRC",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_r <- (sum(CNV_CRC_r))/(ncol(CNV_CRC_r)*nrow(CNV_CRC_r))
      
      CNV_STAD_r <- CNV_r[grep("STAD",rownames(CNV_r)),grep("STAD",colnames(CNV_r))]
      average_CNV_STAD_r <- (sum(CNV_STAD_r))/(ncol(CNV_STAD_r)*nrow(CNV_STAD_r))
      
      CNV_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("NC",colnames(CNV_r))]
      average_CNV_NC_r <- (sum(CNV_NC_r))/(ncol(CNV_NC_r)*nrow(CNV_NC_r))
      
      CNV_CRC_STAD_r <- CNV_r[grep("STAD",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_STAD_r <- (sum(CNV_CRC_STAD_r))/(ncol(CNV_CRC_STAD_r)*nrow(CNV_CRC_STAD_r))
      
      CNV_CRC_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_NC_r <- (sum(CNV_CRC_NC_r))/(ncol(CNV_CRC_NC_r)*nrow(CNV_CRC_NC_r))
      
      CNV_STAD_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("STAD",colnames(CNV_r))]
      average_CNV_STAD_NC_r <- (sum(CNV_STAD_NC_r))/(ncol(CNV_STAD_NC_r)*nrow(CNV_STAD_NC_r))
      
      
      ##### RNA vs. Methylation rectangle
      Expression_CRC_Methy_CRC_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_CRC_r <- (sum(Expression_CRC_Methy_CRC_r))/(ncol(Expression_CRC_Methy_CRC_r)*nrow(Expression_CRC_Methy_CRC_r))
      
      Expression_STAD_Methy_CRC_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_CRC_r <- (sum(Expression_STAD_Methy_CRC_r))/(ncol(Expression_STAD_Methy_CRC_r)*nrow(Expression_STAD_Methy_CRC_r))
      
      Expression_NC_Methy_CRC_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_CRC_r <- (sum(Expression_NC_Methy_CRC_r))/(ncol(Expression_NC_Methy_CRC_r)*nrow(Expression_NC_Methy_CRC_r))
      
      Expression_CRC_Methy_STAD_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_STAD_r <- (sum(Expression_CRC_Methy_STAD_r))/(ncol(Expression_CRC_Methy_STAD_r)*nrow(Expression_CRC_Methy_STAD_r))
      
      Expression_STAD_Methy_STAD_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_STAD_r <- (sum(Expression_STAD_Methy_STAD_r))/(ncol(Expression_STAD_Methy_STAD_r)*nrow(Expression_STAD_Methy_STAD_r))
      
      Expression_NC_Methy_STAD_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_STAD_r <- (sum(Expression_NC_Methy_STAD_r))/(ncol(Expression_NC_Methy_STAD_r)*nrow(Expression_NC_Methy_STAD_r))
      
      Expression_CRC_Methy_NC_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_NC_r <- (sum(Expression_CRC_Methy_NC_r))/(ncol(Expression_CRC_Methy_NC_r)*nrow(Expression_CRC_Methy_NC_r))
      
      Expression_STAD_Methy_NC_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_NC_r <- (sum(Expression_STAD_Methy_NC_r))/(ncol(Expression_STAD_Methy_NC_r)*nrow(Expression_STAD_Methy_NC_r))
      
      Expression_NC_Methy_NC_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_NC_r <- (sum(Expression_NC_Methy_NC_r))/(ncol(Expression_NC_Methy_NC_r)*nrow(Expression_NC_Methy_NC_r))
      ####### RNA vs. CNV rectangle
      Expression_CRC_CNV_CRC_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_CRC_r <- (sum(Expression_CRC_CNV_CRC_r))/(ncol(Expression_CRC_CNV_CRC_r)*nrow(Expression_CRC_CNV_CRC_r))
      
      Expression_STAD_CNV_CRC_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_CRC_r <- (sum(Expression_STAD_CNV_CRC_r))/(ncol(Expression_STAD_CNV_CRC_r)*nrow(Expression_STAD_CNV_CRC_r))
      
      Expression_NC_CNV_CRC_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_CRC_r <- (sum(Expression_NC_CNV_CRC_r))/(ncol(Expression_NC_CNV_CRC_r)*nrow(Expression_NC_CNV_CRC_r))
      
      Expression_CRC_CNV_STAD_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_STAD_r <- (sum(Expression_CRC_CNV_STAD_r))/(ncol(Expression_CRC_CNV_STAD_r)*nrow(Expression_CRC_CNV_STAD_r))
      
      Expression_STAD_CNV_STAD_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_STAD_r <- (sum(Expression_STAD_CNV_STAD_r))/(ncol(Expression_STAD_CNV_STAD_r)*nrow(Expression_STAD_CNV_STAD_r))
      
      Expression_NC_CNV_STAD_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_STAD_r <- (sum(Expression_NC_CNV_STAD_r))/(ncol(Expression_NC_CNV_STAD_r)*nrow(Expression_NC_CNV_STAD_r))
      
      Expression_CRC_CNV_NC_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_NC_r <- (sum(Expression_CRC_CNV_NC_r))/(ncol(Expression_CRC_CNV_NC_r)*nrow(Expression_CRC_CNV_NC_r))
      
      Expression_STAD_CNV_NC_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_NC_r <- (sum(Expression_STAD_CNV_NC_r))/(ncol(Expression_STAD_CNV_NC_r)*nrow(Expression_STAD_CNV_NC_r))
      
      Expression_NC_CNV_NC_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_NC_r <- (sum(Expression_NC_CNV_NC_r))/(ncol(Expression_NC_CNV_NC_r)*nrow(Expression_NC_CNV_NC_r))
      
      ####### Methylation vs CNV rectangle
      Methy_CRC_CNV_CRC_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_CRC_r <- (sum(Methy_CRC_CNV_CRC_r))/(ncol(Methy_CRC_CNV_CRC_r)*nrow(Methy_CRC_CNV_CRC_r))
      
      Methy_STAD_CNV_CRC_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_CRC_r <- (sum(Methy_STAD_CNV_CRC_r))/(ncol(Methy_STAD_CNV_CRC_r)*nrow(Methy_STAD_CNV_CRC_r))
      
      Methy_NC_CNV_CRC_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_CRC_r <- (sum(Methy_NC_CNV_CRC_r))/(ncol(Methy_NC_CNV_CRC_r)*nrow(Methy_NC_CNV_CRC_r))
      
      Methy_CRC_CNV_STAD_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_STAD_r <- (sum(Methy_CRC_CNV_STAD_r))/(ncol(Methy_CRC_CNV_STAD_r)*nrow(Methy_CRC_CNV_STAD_r))
      
      Methy_STAD_CNV_STAD_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_STAD_r <- (sum(Methy_STAD_CNV_STAD_r))/(ncol(Methy_STAD_CNV_STAD_r)*nrow(Methy_STAD_CNV_STAD_r))
      
      Methy_NC_CNV_STAD_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_STAD_r <- (sum(Methy_NC_CNV_STAD_r))/(ncol(Methy_NC_CNV_STAD_r)*nrow(Methy_NC_CNV_STAD_r))
      
      Methy_CRC_CNV_NC_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_NC_r <- (sum(Methy_CRC_CNV_NC_r))/(ncol(Methy_CRC_CNV_NC_r)*nrow(Methy_CRC_CNV_NC_r))
      
      Methy_STAD_CNV_NC_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_NC_r <- (sum(Methy_STAD_CNV_NC_r))/(ncol(Methy_STAD_CNV_NC_r)*nrow(Methy_STAD_CNV_NC_r))
      
      Methy_NC_CNV_NC_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_NC_r <- (sum(Methy_NC_CNV_NC_r))/(ncol(Methy_NC_CNV_NC_r)*nrow(Methy_NC_CNV_NC_r))
      
      average_r <- as.data.frame(matrix(numeric(0),ncol=9,nrow=9))
      
      colnames(average_r) <- c("Expression-CRC","Expression-STAD","Expression-NC","Methylation-CRC","Methylation-STAD","Methylation-NC","CNV-CRC","CNV-STAD","CNV-NC")
      rownames(average_r) <- c("Expression-CRC","Expression-STAD","Expression-NC","Methylation-CRC","Methylation-STAD","Methylation-NC","CNV-CRC","CNV-STAD","CNV-NC")
      
      #Expression triangle
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_CRC_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_CRC_NC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_NC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_STAD_NC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_NC_r
      #Methylation triangle
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_CRC_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_STAD_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_CRC_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_NC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Methy_STAD_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_NC_r
      #copynumber triangle
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-STAD")] <- CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-NC")] <- CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-STAD")] <- CNV_CRC_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_STAD_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-NC")] <- CNV_CRC_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_NC_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-NC")] <- CNV_STAD_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-STAD")] <- CNV_STAD_NC_r
      #Expression vs. Methylation
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_CRC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_CRC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Expression_CRC_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_NC_r
      
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_STAD_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_STAD_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Expression_STAD_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_NC_r
      
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_NC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_NC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-NC")] <- Expression_NC_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_NC_r
      
      #Expression vs. CNV
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-CRC")] <- Expression_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-STAD")] <- Expression_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-NC")] <- Expression_CRC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-CRC")] <- Expression_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-STAD")] <- Expression_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-NC")] <- Expression_STAD_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-CRC")] <- Expression_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-STAD")] <- Expression_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-NC")] <- Expression_NC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_NC_r
      
      #Methylation vs. CNV
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-CRC")] <- Methy_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-STAD")] <- Methy_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-NC")] <- Methy_CRC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-CRC")] <- Methy_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-STAD")] <- Methy_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-NC")] <- Methy_STAD_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-CRC")] <- Methy_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-STAD")] <- Methy_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-NC")] <- Methy_NC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_NC_r
      
      write.csv(average_r,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average-20211026-p.csv")
      #library(corrplot)
      #corrplot(corr = as.matrix(average_r),
      #         method = c("circle"),col = col2(200),
      #         outline = TRUE,tl.col="black",
      #         type="full", 
      #         order="original",
      #         tl.pos = "l",tl.cex=1.5,
      #         cl.cex = 1,cl.lim= c(-0.1,1),cl.length = 12,
      #         number.cex = as.matrix(average_r),
      #         #order = "hclust",
      #         addrect = 3)
    }
    {
      all_r <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired-20211007.csv",header = TRUE,row.names = 1)
      colnames(all_r) <- gsub(".","-",fixed = TRUE,colnames(all_r))
      Expression_r <- all_r[grep("expression",rownames(all_r)),grep("expression",colnames(all_r))]
      Methy_r <- all_r[grep("methylation",rownames(all_r)),grep("methylation",colnames(all_r))]
      CNV_r <- all_r[grep("cnv",rownames(all_r)),grep("cnv",colnames(all_r))]
      Expression_Methy_r <- all_r[grep("expression",rownames(all_r)),grep("methylation",colnames(all_r))]
      Expression_CNV_r <- all_r[grep("expression",rownames(all_r)),grep("cnv",colnames(all_r))]
      Methy_CNV_r <- all_r[grep("methylation",rownames(all_r)),grep("cnv",colnames(all_r))]
      
      ######Expression triangle
      Expression_CRC_r <- Expression_r[grep("CRC",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_r <- (sum(Expression_CRC_r))/(ncol(Expression_CRC_r)*nrow(Expression_CRC_r))
      
      Expression_STAD_r <- Expression_r[grep("STAD",rownames(Expression_r)),grep("STAD",colnames(Expression_r))]
      average_Expression_STAD_r <- (sum(Expression_STAD_r))/(ncol(Expression_STAD_r)*nrow(Expression_STAD_r))
      
      Expression_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("NC",colnames(Expression_r))]
      average_Expression_NC_r <- (sum(Expression_NC_r))/(ncol(Expression_NC_r)*nrow(Expression_NC_r))
      
      Expression_CRC_STAD_r <- Expression_r[grep("STAD",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_STAD_r <- (sum(Expression_CRC_STAD_r))/(ncol(Expression_CRC_STAD_r)*nrow(Expression_CRC_STAD_r))
      
      Expression_CRC_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_NC_r <- (sum(Expression_CRC_NC_r))/(ncol(Expression_CRC_NC_r)*nrow(Expression_CRC_NC_r))
      
      Expression_STAD_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("STAD",colnames(Expression_r))]
      average_Expression_STAD_NC_r <- (sum(Expression_STAD_NC_r))/(ncol(Expression_STAD_NC_r)*nrow(Expression_STAD_NC_r))
      
      ########Methylation triangle
      Methy_CRC_r <- Methy_r[grep("CRC",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_r <- (sum(Methy_CRC_r))/(ncol(Methy_CRC_r)*nrow(Methy_CRC_r))
      
      Methy_STAD_r <- Methy_r[grep("STAD",rownames(Methy_r)),grep("STAD",colnames(Methy_r))]
      average_Methy_STAD_r <- (sum(Methy_STAD_r))/(ncol(Methy_STAD_r)*nrow(Methy_STAD_r))
      
      Methy_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("NC",colnames(Methy_r))]
      average_Methy_NC_r <- (sum(Methy_NC_r))/(ncol(Methy_NC_r)*nrow(Methy_NC_r))
      
      Methy_CRC_STAD_r <- Methy_r[grep("STAD",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_STAD_r <- (sum(Methy_CRC_STAD_r))/(ncol(Methy_CRC_STAD_r)*nrow(Methy_CRC_STAD_r))
      
      Methy_CRC_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_NC_r <- (sum(Methy_CRC_NC_r))/(ncol(Methy_CRC_NC_r)*nrow(Methy_CRC_NC_r))
      
      Methy_STAD_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("STAD",colnames(Methy_r))]
      average_Methy_STAD_NC_r <- (sum(Methy_STAD_NC_r))/(ncol(Methy_STAD_NC_r)*nrow(Methy_STAD_NC_r))
      
      ########CNV triangle
      CNV_CRC_r <- CNV_r[grep("CRC",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_r <- (sum(CNV_CRC_r))/(ncol(CNV_CRC_r)*nrow(CNV_CRC_r))
      
      CNV_STAD_r <- CNV_r[grep("STAD",rownames(CNV_r)),grep("STAD",colnames(CNV_r))]
      average_CNV_STAD_r <- (sum(CNV_STAD_r))/(ncol(CNV_STAD_r)*nrow(CNV_STAD_r))
      
      CNV_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("NC",colnames(CNV_r))]
      average_CNV_NC_r <- (sum(CNV_NC_r))/(ncol(CNV_NC_r)*nrow(CNV_NC_r))
      
      CNV_CRC_STAD_r <- CNV_r[grep("STAD",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_STAD_r <- (sum(CNV_CRC_STAD_r))/(ncol(CNV_CRC_STAD_r)*nrow(CNV_CRC_STAD_r))
      
      CNV_CRC_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_NC_r <- (sum(CNV_CRC_NC_r))/(ncol(CNV_CRC_NC_r)*nrow(CNV_CRC_NC_r))
      
      CNV_STAD_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("STAD",colnames(CNV_r))]
      average_CNV_STAD_NC_r <- (sum(CNV_STAD_NC_r))/(ncol(CNV_STAD_NC_r)*nrow(CNV_STAD_NC_r))
      
      
      ##### RNA vs. Methylation rectangle
      Expression_CRC_Methy_CRC_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_CRC_r <- (sum(Expression_CRC_Methy_CRC_r))/(ncol(Expression_CRC_Methy_CRC_r)*nrow(Expression_CRC_Methy_CRC_r))
      
      Expression_STAD_Methy_CRC_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_CRC_r <- (sum(Expression_STAD_Methy_CRC_r))/(ncol(Expression_STAD_Methy_CRC_r)*nrow(Expression_STAD_Methy_CRC_r))
      
      Expression_NC_Methy_CRC_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_CRC_r <- (sum(Expression_NC_Methy_CRC_r))/(ncol(Expression_NC_Methy_CRC_r)*nrow(Expression_NC_Methy_CRC_r))
      
      Expression_CRC_Methy_STAD_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_STAD_r <- (sum(Expression_CRC_Methy_STAD_r))/(ncol(Expression_CRC_Methy_STAD_r)*nrow(Expression_CRC_Methy_STAD_r))
      
      Expression_STAD_Methy_STAD_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_STAD_r <- (sum(Expression_STAD_Methy_STAD_r))/(ncol(Expression_STAD_Methy_STAD_r)*nrow(Expression_STAD_Methy_STAD_r))
      
      Expression_NC_Methy_STAD_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_STAD_r <- (sum(Expression_NC_Methy_STAD_r))/(ncol(Expression_NC_Methy_STAD_r)*nrow(Expression_NC_Methy_STAD_r))
      
      Expression_CRC_Methy_NC_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_NC_r <- (sum(Expression_CRC_Methy_NC_r))/(ncol(Expression_CRC_Methy_NC_r)*nrow(Expression_CRC_Methy_NC_r))
      
      Expression_STAD_Methy_NC_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_NC_r <- (sum(Expression_STAD_Methy_NC_r))/(ncol(Expression_STAD_Methy_NC_r)*nrow(Expression_STAD_Methy_NC_r))
      
      Expression_NC_Methy_NC_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_NC_r <- (sum(Expression_NC_Methy_NC_r))/(ncol(Expression_NC_Methy_NC_r)*nrow(Expression_NC_Methy_NC_r))
      ####### RNA vs. CNV rectangle
      Expression_CRC_CNV_CRC_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_CRC_r <- (sum(Expression_CRC_CNV_CRC_r))/(ncol(Expression_CRC_CNV_CRC_r)*nrow(Expression_CRC_CNV_CRC_r))
      
      Expression_STAD_CNV_CRC_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_CRC_r <- (sum(Expression_STAD_CNV_CRC_r))/(ncol(Expression_STAD_CNV_CRC_r)*nrow(Expression_STAD_CNV_CRC_r))
      
      Expression_NC_CNV_CRC_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_CRC_r <- (sum(Expression_NC_CNV_CRC_r))/(ncol(Expression_NC_CNV_CRC_r)*nrow(Expression_NC_CNV_CRC_r))
      
      Expression_CRC_CNV_STAD_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_STAD_r <- (sum(Expression_CRC_CNV_STAD_r))/(ncol(Expression_CRC_CNV_STAD_r)*nrow(Expression_CRC_CNV_STAD_r))
      
      Expression_STAD_CNV_STAD_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_STAD_r <- (sum(Expression_STAD_CNV_STAD_r))/(ncol(Expression_STAD_CNV_STAD_r)*nrow(Expression_STAD_CNV_STAD_r))
      
      Expression_NC_CNV_STAD_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_STAD_r <- (sum(Expression_NC_CNV_STAD_r))/(ncol(Expression_NC_CNV_STAD_r)*nrow(Expression_NC_CNV_STAD_r))
      
      Expression_CRC_CNV_NC_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_NC_r <- (sum(Expression_CRC_CNV_NC_r))/(ncol(Expression_CRC_CNV_NC_r)*nrow(Expression_CRC_CNV_NC_r))
      
      Expression_STAD_CNV_NC_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_NC_r <- (sum(Expression_STAD_CNV_NC_r))/(ncol(Expression_STAD_CNV_NC_r)*nrow(Expression_STAD_CNV_NC_r))
      
      Expression_NC_CNV_NC_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_NC_r <- (sum(Expression_NC_CNV_NC_r))/(ncol(Expression_NC_CNV_NC_r)*nrow(Expression_NC_CNV_NC_r))
      
      ####### Methylation vs CNV rectangle
      Methy_CRC_CNV_CRC_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_CRC_r <- (sum(Methy_CRC_CNV_CRC_r))/(ncol(Methy_CRC_CNV_CRC_r)*nrow(Methy_CRC_CNV_CRC_r))
      
      Methy_STAD_CNV_CRC_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_CRC_r <- (sum(Methy_STAD_CNV_CRC_r))/(ncol(Methy_STAD_CNV_CRC_r)*nrow(Methy_STAD_CNV_CRC_r))
      
      Methy_NC_CNV_CRC_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_CRC_r <- (sum(Methy_NC_CNV_CRC_r))/(ncol(Methy_NC_CNV_CRC_r)*nrow(Methy_NC_CNV_CRC_r))
      
      Methy_CRC_CNV_STAD_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_STAD_r <- (sum(Methy_CRC_CNV_STAD_r))/(ncol(Methy_CRC_CNV_STAD_r)*nrow(Methy_CRC_CNV_STAD_r))
      
      Methy_STAD_CNV_STAD_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_STAD_r <- (sum(Methy_STAD_CNV_STAD_r))/(ncol(Methy_STAD_CNV_STAD_r)*nrow(Methy_STAD_CNV_STAD_r))
      
      Methy_NC_CNV_STAD_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_STAD_r <- (sum(Methy_NC_CNV_STAD_r))/(ncol(Methy_NC_CNV_STAD_r)*nrow(Methy_NC_CNV_STAD_r))
      
      Methy_CRC_CNV_NC_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_NC_r <- (sum(Methy_CRC_CNV_NC_r))/(ncol(Methy_CRC_CNV_NC_r)*nrow(Methy_CRC_CNV_NC_r))
      
      Methy_STAD_CNV_NC_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_NC_r <- (sum(Methy_STAD_CNV_NC_r))/(ncol(Methy_STAD_CNV_NC_r)*nrow(Methy_STAD_CNV_NC_r))
      
      Methy_NC_CNV_NC_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_NC_r <- (sum(Methy_NC_CNV_NC_r))/(ncol(Methy_NC_CNV_NC_r)*nrow(Methy_NC_CNV_NC_r))
      
      average_r <- as.data.frame(matrix(numeric(0),ncol=9,nrow=9))
      
      colnames(average_r) <- c("Expression-CRC","Expression-STAD","Expression-NC","Methylation-CRC","Methylation-STAD","Methylation-NC","CNV-CRC","CNV-STAD","CNV-NC")
      rownames(average_r) <- c("Expression-CRC","Expression-STAD","Expression-NC","Methylation-CRC","Methylation-STAD","Methylation-NC","CNV-CRC","CNV-STAD","CNV-NC")
      
      #Expression triangle
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_CRC_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_CRC_NC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_NC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_STAD_NC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_NC_r
      #Methylation triangle
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_CRC_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_STAD_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_CRC_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_NC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Methy_STAD_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_NC_r
      #copynumber triangle
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-STAD")] <- CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-NC")] <- CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-STAD")] <- CNV_CRC_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_STAD_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-NC")] <- CNV_CRC_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_NC_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-NC")] <- CNV_STAD_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-STAD")] <- CNV_STAD_NC_r
      #Expression vs. Methylation
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_CRC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_CRC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Expression_CRC_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_NC_r
      
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_STAD_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_STAD_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Expression_STAD_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_NC_r
      
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_NC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_NC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-NC")] <- Expression_NC_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_NC_r
      
      #Expression vs. CNV
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-CRC")] <- Expression_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-STAD")] <- Expression_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-NC")] <- Expression_CRC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-CRC")] <- Expression_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-STAD")] <- Expression_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-NC")] <- Expression_STAD_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-CRC")] <- Expression_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-STAD")] <- Expression_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-NC")] <- Expression_NC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_NC_r
      
      #Methylation vs. CNV
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-CRC")] <- Methy_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-STAD")] <- Methy_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-NC")] <- Methy_CRC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-CRC")] <- Methy_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-STAD")] <- Methy_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-NC")] <- Methy_STAD_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-CRC")] <- Methy_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-STAD")] <- Methy_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-NC")] <- Methy_NC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_NC_r
      
      write.csv(average_r,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average-20211026.csv")
      #library(corrplot)
      #corrplot(corr = as.matrix(average_r),
      #         method = c("circle"),col = col2(200),
      #         outline = TRUE,tl.col="black",
      #         type="full", 
      #         order="original",
      #         tl.pos = "l",tl.cex=1.5,
      #         cl.cex = 1,cl.lim= c(-0.1,1),cl.length = 12,
      #         number.cex = as.matrix(average_r),
      #         #order = "hclust",
      #         addrect = 3)
    }
    
    r <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average-20211026.csv",header = TRUE,row.names = 1)
    p <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average-20211026-p.csv",header = TRUE, row.names = 1)
    p[is.na(p)] <- 0
    r$omics <- rownames(r)
    r_forplot <- melt(r)
    r_forplot$variable <- gsub(".","-",fixed=TRUE, r_forplot$variable)
    #paste0(r_forplot$variable,r_forplot$omics)
    colnames(r_forplot) <- c("omics1","omics2","pearson")
    
    p$omics <- rownames(p)
    p_forplot <- melt(p)
    p_forplot$variable <- gsub(".","-",fixed=TRUE, p_forplot$variable)
    #paste0(r_forplot$variable,r_forplot$omics)
    colnames(p_forplot) <- c("omics3","omics4","pvalue")
    corplot <- cbind(r_forplot,p_forplot)
    corplot <- corplot[,-which(colnames(corplot) == "omics3")]
    corplot <- corplot[,-which(colnames(corplot) == "omics4")]
    
    corplot$omics1 <- factor(corplot$omics1,c("Expression-STAD","Expression-CRC","Expression-NC",
                                              "Methylation-STAD","Methylation-CRC","Methylation-NC",
                                              "CNV-STAD","CNV-CRC","CNV-NC"))
    corplot$omics2 <- factor(corplot$omics2,rev(c("Expression-STAD","Expression-CRC","Expression-NC",
                                              "Methylation-STAD","Methylation-CRC","Methylation-NC",
                                              "CNV-STAD","CNV-CRC","CNV-NC")))
    ggplot(corplot,aes(x=omics1,y=omics2))+
      geom_point(aes(size=as.numeric(pearson),color=as.numeric(pvalue)))+
      geom_vline(xintercept = c(3.5,6.5))+
      geom_vline(xintercept = c(1.5,2.5,4.5,5.5,7.5,8.5),colour = "grey",size = 0.1)+
      geom_hline(yintercept = c(3.5,6.5))+
      geom_hline(yintercept = c(1.5,2.5,4.5,5.5,7.5,8.5),colour = "grey",size = 0.1)+
      scale_colour_gradient2(low="#B6202E",high="#549EC9",mid="white", midpoint = 0.5)+
      #geom_vline(xintercept = c(6.5,10.5,16.5),linetype = "dashed",color="black")+
      #geom_vline(xintercept = 25.5,linetype = "dashed",color="grey")+
      labs(color=expression(PValue),
           size="Pearson's correlation",
           x="")+
      theme_bw()+
      theme(panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
            axis.text.x = element_text(size=rel(1.3),face="bold",colour = "black",angle = 45,vjust=1,hjust=1),
            axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
            axis.title.y = element_blank())
  }
  
  #plasma vs PBMC vs tumor(not in multi-omics paper)
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

#normalized gene coverage between RNA and DNA: STAD GAPDH
{
  {
    {
      library(ggbreak) 
      library(patchwork)
      target_gene <- "GAPDH"
      setwd(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Coverage/STAD-PKU-13_normalized/",target_gene))
      DNA <- read.csv("STAD-PKU-13-wgs_sorted_cpm.txt",sep = "\t",header = FALSE)
      methy <- read.csv("STAD-PKU-13-me_sorted_cpm.txt",sep = "\t",header = FALSE)
      #RNA <- read.csv("STAD-PKU-13-pico-intron-spanning_sorted_cpm.txt",sep = "\t",header = FALSE)
      RNA <- read.csv("STAD-PKU-13-genome_rmdup_sorted_cpm.txt",sep = "\t",header = FALSE)
      
      if(file.info("STAD-PKU-13-qia_sorted_cpm.txt")$size == 0) {
        qia <- as.data.frame(matrix(numeric(0),ncol=3))
      } else {
        qia <- read.csv("STAD-PKU-13-qia_sorted_cpm.txt",sep = "\t",header = FALSE)
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
      #plot <- plot[(plot$Position>7668000)&(plot$Position<7687550),]
      
      RNA_upperlimit <- ceiling(max(plot[grep("Total cfRNA",plot$type),]$Depth)/0.1)*0.1
      WGS_meanCoverage <- ceiling(mean(plot[grep("cfDNA WGS",plot$type),]$Depth)/0.1)*0.1
      #WGS_meanCoverage <- 20
      
      p <- ggplot(plot,aes(x=plot$Position,y=log(plot$Depth+1,10),fill=plot$type))+
        geom_bar(stat="identity",position = 'dodge')+
        #geom_bar(aes(x=plot$Position,y=log10(plot$layer2+1),fill=plot$type),stat="identity",position = 'dodge')+
        geom_bar(aes(x=plot$Position,y=log10(plot$layer4+1),fill=plot$type),stat="identity",position = 'dodge')+
        geom_bar(aes(x=plot$Position,y=log10(plot$layer3+1),fill=plot$type),stat="identity",position = 'dodge')+
        scale_x_continuous(expand = c(0,0))+
        scale_y_continuous(limits=c(0, log(1+RNA_upperlimit+RNA_upperlimit*0.8,10)),
                           breaks = c(0,log(WGS_meanCoverage+1,10),log(RNA_upperlimit+1,10)),
                           labels=c("",WGS_meanCoverage,RNA_upperlimit),expand = c(0,0))+
        #scale_x_discrete(expand=c(0.2, 0.2))+
        #scale_fill_jco(alpha = 0.8)+
        #scale_fill_manual(values=c(alpha("#5CACEE","#FFCC11",alpha = 0.45),"#666666","#CD0000"))+
        scale_fill_manual(values=c(alpha("#5CACEE",alpha = 1),alpha("#EEB4B4",alpha = 0.45),alpha("#666666",alpha = 0.1),alpha("#CD0000",alpha = 0.8)))+
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
      
      p <- p+scale_y_break(c(log(WGS_meanCoverage*10,10), log(1+RNA_upperlimit-RNA_upperlimit*0.8,10)),
                           scale = 0.7,ticklabels = NULL)+
        scale_y_continuous(limits=c(0, log(1+RNA_upperlimit+RNA_upperlimit*0.8,10)),
                           breaks = c(0,log(WGS_meanCoverage+1,10),log(RNA_upperlimit-100+1,10)),
                           labels=c("",WGS_meanCoverage,RNA_upperlimit-100),expand = c(0,0))
      
      ggsave(paste0(target_gene,"_genome_rmdup.pdf"),p,width=18,height=6)
    }
    
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(trackViewer)
    #GAPDH: chr12	6533926	6538374	ENSG00000111640.14	0	+
    #APC: chr5	112707497	112846239	ENSG00000134982.16	0	+
    #KRAS: chr12	25204788	25250936	ENSG00000133703.11	0	-
    #TP53: chr17	7661778	7687550	ENSG00000141510.16	0	-
    gr <- GRanges("chr17", IRanges(7668000, 7687550), strand="+")
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
