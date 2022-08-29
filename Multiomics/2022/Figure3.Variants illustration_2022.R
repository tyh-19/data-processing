library(rsvg)
rsvg_pdf("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/15.examples/expression_20220209/MIR629.svg", file = "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/15.examples/expression_20220209/MIR629.pdf")

setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/15.examples/qPCR/")


#longRNA (2 group comparison)
{
#read in qPCR results
  {
  qPCR_result <- read.csv(file = "alllong_gapdhnorm_1124_plot3_reshaped.csv",header = TRUE)
  qPCR_result$Group <- factor(qPCR_result$Group,levels=c("STAD","Gastritis","HD"))
  qPCR_result$Candidates <- factor(qPCR_result$Candidates,
                                   levels=c("ADAMTS6-2",
                                            "ADAMTS6-5",
                                            "CORO1A-1",
                                            "Lnc-COL6A3-1",
                                            "GClnc1",
                                            "MT-ATP6",
                                            "MT-TP",
                                            "S100A4",
                                            "SERPING1-1",
                                            "SSX2IP-1",
                                            "SSX2IP-2",
                                            "SSX2IP-8",
                                            "TRAPPC1"))
  qPCR_result$Candidates_Group <- factor(qPCR_result$Candidates_Group,
                                         levels=c("ADAMTS6-2_STAD","ADAMTS6-2_HD",
                                                  "ADAMTS6-5_STAD","ADAMTS6-5_HD",
                                                  "CORO1A-1_STAD","CORO1A-1_HD",
                                                  "Lnc-COL6A3-1_STAD","Lnc-COL6A3-1_HD",
                                                  "GClnc1_STAD","GClnc1_HD",
                                                  "MT-ATP6_STAD","MT-ATP6_HD",
                                                  "MT-TP_STAD","MT-TP_HD",
                                                  "S100A4_STAD","S100A4_HD",
                                                  "SERPING1-1_STAD","SERPING1-1_HD",
                                                  "SSX2IP-1_STAD","SSX2IP-1_HD",
                                                  "SSX2IP-2_STAD","SSX2IP-2_HD",
                                                  "SSX2IP-8_STAD","SSX2IP-8_HD",
                                                  "TRAPPC1_STAD","TRAPPC1_HD"))
  }
  
#plot
{
#main_plot
{
p <- ggplot(qPCR_result,aes(x=Candidates_Group,y=neg_deltaCt,colour=Group, width = 5),fullrange = FALSE)+
  geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
  geom_point(size = 2.5, position = position_jitterdodge(dodge.width=1,jitter.width = 0.1))+
  #geom_point(aes(shape=Hospital),size = 3, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
  theme_bw()+
  theme(plot.margin = unit(c(1, 1, 4, 1),"cm"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=30),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=30),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=30,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=36),
        axis.title.y = element_text(face="bold",color="black", size=36))+
  coord_cartesian(clip="off")+#key for annotate candidate 
  #https://stackoverflow.com/questions/64654407/annotate-ggplot2-across-both-axis-text-keeps-changing-position
  labs(x="",y="-delta Ct",title="", face="bold")+scale_color_manual(values = c("STAD"="#9D1309","HD"="#3A5894","Gastritis"="#556B2F"))
}

#add line
{
  group_num=length(unique(qPCR_result$Group))
  #https://stackoverflow.com/questions/44317502/update-a-ggplot-using-a-for-loop-r
  
  i=1
  while(i+1<=length(unique(qPCR_result$Candidates))){
    p <- p+geom_vline(aes_(xintercept = 0.5+i*group_num),linetype="dashed",size = 0.8)  
    ## group_num controls vline's position
    i=i+1
  }
  }

#get y limits for FC label
{
  a <- as.list(ggplot_build(p)$layout)
  b <- unlist(a$panel_params)
  y_bottom <- b$y.range1
  y_up <- b$y.range2
}

#wilcox test by group, label FC and x.axis
##test between each 2 columns, positive:1st column,negative:2n columm. When not sure, check groups.
{
  group_num=length(unique(qPCR_result$Group))  
  
  i=1
  j=1
  groups <- as.character(levels(qPCR_result$Candidates_Group))
  xlabel <- as.character(levels(qPCR_result$Candidates))
  while(i<=length(unique(qPCR_result$Candidates_Group))){
    FC <- mean(subset(qPCR_result,Candidates_Group==groups[i])$neg_deltaCt,na.rm = TRUE)-mean(subset(qPCR_result,Candidates_Group==groups[i+1])$neg_deltaCt,na.rm = TRUE)
    #https://www.jianshu.com/p/aab34be5f983
    k=1
    l <- list()
    while(k<group_num){
      l_tmp <- list(c(groups[i],groups[i+k]))
      if(mean(subset(qPCR_result,Candidates_Group==l_tmp[[1]][1])$neg_deltaCt,na.rm = TRUE)=="NaN" ||
         mean(subset(qPCR_result,Candidates_Group==l_tmp[[1]][2])$neg_deltaCt,na.rm = TRUE)=="NaN"){
        k=k+1
      }
      else {
        l <- c(l,l_tmp)
        k=k+1
      }
    }
    if(FC >= 0 || FC=="NaN"){
      p <- p+stat_compare_means(comparison=l,
                                label = "p.signif",
                                hide.ns = TRUE,
                                size=8,
                                vjust=0.5,
                                method="wilcox.test",
                                method.args = list(alternative = "greater"))+
        #annotate("text",size=5,x=i+0.5,y=y_up-7,
        #         label=paste0("-△△Ct:",round(FC,2),"\nWilcox\nGreater"),family = "Arial", fontface = "bold")+
        #解释annotate: https://blog.csdn.net/g_r_c/article/details/19673625
        annotate("text",size=8,x=i+group_num/2-1+0.5,y=-Inf,
                 label=xlabel[j],family = "Arial", fontface = "bold",angle=45,hjust=1,vjust=1)
      #geom_text(aes_(x=i+0.5,y=10,label=paste0("FC:",round(FC,3))))
    } else {
      p <- p+stat_compare_means(comparison=l,
                                label = "p.signif",
                                hide.ns = TRUE,
                                size=8,
                                vjust=0.5,
                                method="wilcox.test",
                                method.args = list(alternative = "less"))+
        #annotate("text",size=5,x=i+0.5,y=y_up-1,
        #         label=paste0("-△△Ct:",round(FC,2),"\nWilcox\nLess"),family = "Arial", fontface = "bold")+
        annotate("text",size=8,x=i+group_num/2-1+0.5,y=-Inf,
                 label=xlabel[j],family = "Arial", fontface = "bold",angle=45,hjust=1,vjust=1)
    }
    i=i+group_num ## group_num controls gene name's position
    j=j+1
  }
}
  p
}
}

#miRNA (3 group comparison, with NA group)
{
#read in qPCR results
{
  qPCR_result <- read.csv(file = "allmiRNA_mir16norm_1112_plot_reshaped.csv",header = TRUE)
  qPCR_result$Group <- factor(qPCR_result$Group,levels=c("STAD","Gastritis","HD"))
  qPCR_result <- qPCR_result[grep("^MIR",qPCR_result$Candidates, invert = TRUE),]
  #qPCR_result <- qPCR_result[rowSums(is.na(qPCR_result[,1:2])) == 0,]
  qPCR_result$Candidates <- factor(qPCR_result$Candidates,
                                         levels=c(#"MIR486",
                                                  #"MIR92B-1",
                                                  #"MIR92B-2",
                                                  #"MIR99B-2",
                                                  "MT-TD",
                                                  "MT-TG",
                                                  "MT-TL2",
                                                  "RN7SL1-3",
                                                  "RN7SL2-1",
                                                  "CORO1A-2F",
                                                  "CORO1A-4",
                                                  "CORO1A-6",
                                                  "SSX2IP-1F"))
  qPCR_result$Candidates_Group <- factor(qPCR_result$Candidates_Group,
                                         levels=c(#"MIR486_STAD","MIR486_Gastritis","MIR486_HD",
                                                  #"MIR92B-1_STAD","MIR92B-1_Gastritis","MIR92B-1_HD",
                                                  #"MIR92B-2_STAD","MIR92B-2_Gastritis","MIR92B-2_HD",
                                                  #"MIR99B-2_STAD","MIR99B-2_Gastritis","MIR99B-2_HD",
                                                  "MT-TD_STAD","MT-TD_Gastritis","MT-TD_HD",
                                                  "MT-TG_STAD","MT-TG_Gastritis","MT-TG_HD",
                                                  "MT-TL2_STAD","MT-TL2_Gastritis","MT-TL2_HD",
                                                  "RN7SL1-3_STAD","RN7SL1-3_Gastritis","RN7SL1-3_HD",
                                                  "RN7SL2-1_STAD","RN7SL2-1_Gastritis","RN7SL2-1_HD",
                                                  "CORO1A-2F_STAD","CORO1A-2F_Gastritis","CORO1A-2F_HD",
                                                  "CORO1A-4_STAD","CORO1A-4_Gastritis","CORO1A-4_HD",
                                                  "CORO1A-6_STAD","CORO1A-6_Gastritis","CORO1A-6_HD",
                                                  "SSX2IP-1F_STAD","SSX2IP-1F_Gastritis","SSX2IP-1F_HD"))
  
}

#plot
{
  #main_plot
  {
    p <- ggplot(qPCR_result,aes(x=Candidates_Group,y=neg_deltaCt,colour=Group, width = 5),fullrange = FALSE)+
      geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
      geom_point(size = 2, position = position_jitterdodge(dodge.width=1,jitter.width = 0.1))+
      #geom_point(aes(shape=Hospital),size = 3, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
      theme_bw()+
      theme(plot.margin = unit(c(1, 1, 4, 1),"cm"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            axis.line = element_line(size=1, colour = "black"),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=30),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=30),
            plot.title = element_text(hjust = 0.5,size=36,face="bold"),
            axis.text.x = element_blank(),
            #axis.text.x = element_text(face="bold", color="black", size=30,angle = 45,hjust = 1),
            axis.text.y = element_text(face="bold",  color="black", size=24),
            axis.title.x = element_text(face="bold", color="black", size=36),
            axis.title.y = element_text(face="bold",color="black", size=36))+
      coord_cartesian(clip="off")+#key for annotate candidate 
      #https://stackoverflow.com/questions/64654407/annotate-ggplot2-across-both-axis-text-keeps-changing-position
      labs(x="",y="-delta Ct",title="", face="bold")+scale_color_manual(values = c("STAD"="#9D1309","HD"="#3A5894","Gastritis"="#556B2F"))
  }
  
  #add line
  {
    group_num=length(unique(qPCR_result$Group))
    #https://stackoverflow.com/questions/44317502/update-a-ggplot-using-a-for-loop-r
    
    i=1
    while(i+1<=length(unique(qPCR_result$Candidates))){
      p <- p+geom_vline(aes_(xintercept = 0.5+i*group_num),linetype="dashed",size = 0.8)  
      ## group_num controls vline's position
      i=i+1
    }
  }
  
  #get y limits for FC label
  {
    a <- as.list(ggplot_build(p)$layout)
    b <- unlist(a$panel_params)
    y_bottom <- b$y.range1
    y_up <- b$y.range2
  }
  
  #wilcox test by group, label FC and x.axis
  ##test between each 2 columns, positive:1st column,negative:2n columm. When not sure, check groups.
  ##set FC == "Two.sided" means all test are two.sided
  {
    #FC <- "Two.sided"
    group_num=length(unique(qPCR_result$Group))  
    
    groups <- as.character(levels(qPCR_result$Candidates_Group))
    xlabel <- as.character(levels(qPCR_result$Candidates))
    k=1
    while(k<group_num){
      i=1
      j=1
      while(i<=length(unique(qPCR_result$Candidates_Group))){
      FC <- mean(subset(qPCR_result,Candidates_Group==groups[i])$neg_deltaCt,na.rm = TRUE)-mean(subset(qPCR_result,Candidates_Group==groups[i+k])$neg_deltaCt,na.rm = TRUE)
      #https://www.jianshu.com/p/aab34be5f983

      l <- list()
        l_tmp <- list(c(groups[i],groups[i+k]))
        if(mean(subset(qPCR_result,Candidates_Group==l_tmp[[1]][1])$neg_deltaCt,na.rm = TRUE)=="NaN" ||
           mean(subset(qPCR_result,Candidates_Group==l_tmp[[1]][2])$neg_deltaCt,na.rm = TRUE)=="NaN"){
          print(paste0("Cannot compare ",paste(l_tmp[[1]][1],l_tmp[[1]][2], sep = " vs. ")))
        }
        else {
        l <- c(l,l_tmp)
        }
      
      if(FC == "Two.sided"){
        p <- p+stat_compare_means(comparison=l,
                                  label = "p.signif",
                                  size=5,
                                  method="wilcox.test",
                                  method.args = list(alternative = "two.sided"))+
          annotate("text",size=8,x=i+group_num/2-1+0.5,y=-Inf,
                   label=xlabel[j],family = "Arial", fontface = "bold",angle=45,hjust=1,vjust=1)
        #geom_text(aes_(x=i+0.5,y=10,label=paste0("FC:",round(FC,3))))
      } else if(FC >= 0 || FC=="NaN"){
        p <- p+stat_compare_means(comparison=l,
                                  vjust=0.5,
                                  hide.ns = TRUE,
                                  label.y = y_up+3*k-5,
                                  label = "p.signif",
                                  size=8,
                                  #color = "red",
                                  method="wilcox.test",
                                  method.args = list(alternative = "greater"))+
          #annotate("text",size=5,x=i+0.5,y=y_up-7,
          #         label=paste0("-△△Ct:",round(FC,2),"\nWilcox\nGreater"),family = "Arial", fontface = "bold")+
          #解释annotate: https://blog.csdn.net/g_r_c/article/details/19673625
          annotate("text",size=8,x=i+group_num/2-1+0.5,y=-Inf,
                   label=xlabel[j],family = "Arial", fontface = "bold",angle=45,hjust=1,vjust=1)
        #geom_text(aes_(x=i+0.5,y=10,label=paste0("FC:",round(FC,3))))
      } else {
        p <- p+stat_compare_means(comparison=l,
                                  vjust=0.5,
                                  hide.ns = TRUE,
                                  label.y = y_up+3*k-5,
                                  label = "p.signif",
                                  size=8,
                                  #color = "red",
                                  method="wilcox.test",
                                  method.args = list(alternative = "less"))+
          #annotate("text",size=5,x=i+0.5,y=y_up-1,
          #         label=paste0("-△△Ct:",round(FC,2),"\nWilcox\nLess"),family = "Arial", fontface = "bold")+
          annotate("text",size=8,x=i+group_num/2-1+0.5,y=-Inf,
                   label=xlabel[j],family = "Arial", fontface = "bold",angle=45,hjust=1,vjust=1)
      }
      i=i+group_num ## group_num controls gene name's position
      j=j+1
      }
    p
    k=k+1
    }  
  }
  p
}
}

#heatmap for Expression
{
  Differential <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/Expression/output/Expression_GIvsNC_edger_exact.txt",header = TRUE, row.names = 1, sep = "\t")
  Differential_TPM <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/Expression/matrix/Expression_for_GIvsNC_TPM.txt",header = TRUE, row.names = 1, sep = "\t")
  Differential$Gene <- rownames(Differential)
  Differential_TPM$Gene <- rownames(Differential_TPM)
  
  forheatmap <- left_join(Differential,Differential_TPM,by=c("Gene"="Gene"))
  rownames(forheatmap) <- forheatmap$Gene
  forheatmap <- forheatmap[,-which(colnames(forheatmap)=="Gene")]
  forheatmap <- forheatmap[grep("ENSG",rownames(forheatmap)),]
  forheatmap <- forheatmap[,-which(colnames(forheatmap)=="CRC.PKU.5.pico")]
  forheatmap <- forheatmap[,-which(colnames(forheatmap)=="NC.PKU.mix17.pico")]
  forheatmap <- forheatmap[,-which(colnames(forheatmap)=="STAD.PKU.4.pico")]
  #forheatmap <- forheatmap[,-which(colnames(forheatmap)=="CRC.PKU.34.pico")]
  
  gini_Expression_STAD <- gini(t(forheatmap[,grep("STAD",colnames(forheatmap))]))
  gini_Expression_STAD <- as.data.frame(gini_Expression_STAD)
  colnames(gini_Expression_STAD) <- c("gini_index_STAD")
  gini_Expression_CRC <- gini(t(forheatmap[,grep("CRC",colnames(forheatmap))]))
  gini_Expression_CRC <- as.data.frame(gini_Expression_CRC)
  colnames(gini_Expression_CRC) <- c("gini_index_CRC")
  gini_Expression_GI <- gini(t(forheatmap[,grep("STAD|CRC",colnames(forheatmap))]))
  gini_Expression_GI <- as.data.frame(gini_Expression_GI)
  colnames(gini_Expression_GI) <- c("gini_index_GI")
  gini_Expression_NC <- gini(t(forheatmap[,grep("NC",colnames(forheatmap))]))
  gini_Expression_NC <- as.data.frame(gini_Expression_NC)
  colnames(gini_Expression_NC) <- c("gini_index_NC")
  gini_Expression <- cbind(gini_Expression_CRC,gini_Expression_STAD,gini_Expression_GI,gini_Expression_NC)
  
  forheatmap <- cbind(gini_Expression, forheatmap)
  
  heatmap <- forheatmap[(forheatmap$padj < 0.05) & (sqrt(forheatmap$log2FoldChange*forheatmap$log2FoldChange) > 1.5) 
                        #& (forheatmap$gini_index_CRC < 0.5) & (forheatmap$gini_index_STAD < 0.5)
                        & (forheatmap$gini_index_GI < 0.6) & (forheatmap$gini_index_NC < 0.6),]
  heatmap <- heatmap[order(heatmap$log2FoldChange),]
  rownames(heatmap) <- as.character(lapply(strsplit(rownames(heatmap),"\\."), function(x) x[1]))
  
  row_annotation <- data.frame(row.names = rownames(heatmap),log2FoldChange = heatmap$log2FoldChange)
  row_annotation$Trend <- ""
  row_annotation[row_annotation$log2FoldChange > 0,]$Trend <- "Up regulated in cancer"
  row_annotation[row_annotation$log2FoldChange < 0,]$Trend <- "Down regulated in cancer"
  row_annotation <- data.frame(row.names = rownames(row_annotation),Trend = row_annotation$Trend)
  row_annotation$Trend <- factor(row_annotation$Trend, levels = c("Down regulated in cancer","Up regulated in cancer"), ordered = TRUE)
  
  col_annotation <- data.frame(row.names = colnames(heatmap),CancerType=as.character(lapply(strsplit(colnames(heatmap),"\\."), function(x) x[1])))
  col_annotation$CancerType <- gsub("NC","HD",col_annotation$CancerType)
  col_annotation$CancerType <- factor(col_annotation$CancerType, levels = c("CRC","STAD","HD"), ordered = TRUE)
  col_annotation <- na.omit(col_annotation)
  heatmap_CRC <- heatmap[,grep("CRC",colnames(heatmap))]
  heatmap_STAD <- heatmap[,grep("STAD",colnames(heatmap))]
  heatmap_NC <- heatmap[,grep("NC",colnames(heatmap))]
  heatmap <- cbind(heatmap_CRC,heatmap_STAD,heatmap_NC)
  
  #library(biomaRt)
  #mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  
  value <- as.character(lapply(strsplit(rownames(heatmap),"\\."),function(x) x[1]))
  
  gene_names <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values=value, mart= mart,useCache = FALSE)
  
  gene_names[which(gene_names$ensembl_gene_id=="ENSG00000223911"),which(colnames(gene_names)=="hgnc_symbol")] <- "ENSG00000223911"
  gene_names[which(gene_names$ensembl_gene_id=="ENSG00000235298"),which(colnames(gene_names)=="hgnc_symbol")] <- "ENSG00000235298"
  
  heatmap$ensembl_gene_id <- as.character(lapply(strsplit(rownames(heatmap),"\\."),function(x) x[1]))
  heatmap <- left_join(heatmap,gene_names,by=c("ensembl_gene_id"="ensembl_gene_id"))
  rownames(heatmap) <- heatmap$hgnc_symbol
  heatmap <- heatmap[,-which(colnames(heatmap)=="ensembl_gene_id")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="hgnc_symbol")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="gini_index_NC")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="gini_index_CRC")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="gini_index_STAD")]
  
  row_annotation$ensembl_gene_id <- rownames(row_annotation)
  row_annotation <- left_join(row_annotation,gene_names,by=c("ensembl_gene_id"="ensembl_gene_id"))
  row_annotation <- data.frame(row.names = row_annotation$hgnc_symbol,Trend = row_annotation$Trend)
  row_annotation$Trend <- factor(row_annotation$Trend, levels = c("Down regulated in cancer","Up regulated in cancer"), ordered = TRUE)
  ann_color <- list(CancerType=c("CRC"="#EE7621","STAD"="red","HD"="blue"),Trend=c("Down regulated in cancer"="#D6D6D6","Up regulated in cancer"="#EB5E66"))
  bk = unique(c(seq(-3,3, length=100)))
  ComplexHeatmap::pheatmap(
    as.matrix(heatmap),
    breaks = bk,
    annotation_col = col_annotation,
    annotation_row = row_annotation,
    annotation_colors=ann_color,
    gaps_row=c(14),
    gaps_col=c(75),
    scale = "row",
    cluster_cols = FALSE,cluster_rows = FALSE,
    show_colnames= FALSE, 
    show_rownames=TRUE,
    #cluster_cols = FALSE,
    colorRampPalette(c("black","#EEE9E9","red"))(100),
    fontsize_row = 10,
    name = "TPM")
  
  }

#miRNA QC
{
Small_RNA_QC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/02.Quality control/Plasma_smallRNA/Small_RNA.csv",header = TRUE,row.names = 1)

Small_RNA_QC$Pass <- ""
Small_RNA_QC[(Small_RNA_QC$miRNA<10^6) | (Small_RNA_QC$miRNA.genome.ratio < 0.5),]$Pass <- "red"
Small_RNA_QC$ID <- rownames(Small_RNA_QC)
Small_RNA_QC$ID <- gsub(".","-",fixed=TRUE,Small_RNA_QC$ID)
Small_RNA_QC$ID <- gsub("_1","",fixed=TRUE,Small_RNA_QC$ID)
Small_RNA_QC$ID <- gsub("-qia","",fixed=TRUE,Small_RNA_QC$ID)
Small_RNA_QC$ID <- gsub("CRC","CRC-PKU",fixed=TRUE,Small_RNA_QC$ID)
Small_RNA_QC[(Small_RNA_QC$miRNA>=10^6) & (Small_RNA_QC$miRNA.genome.ratio >= 0.5),]$ID <- ""

library(ggrepel)
ggplot(Small_RNA_QC)+geom_point(aes(x=log10(miRNA), y=miRNA.genome.ratio, color = Pass))+
  scale_color_manual(values = c("black","red"))+
  geom_text_repel(aes(x=log10(miRNA), y=miRNA.genome.ratio,label=ID))+
  geom_vline(xintercept = log10(10^6), linetype="dotdash",color = "grey")+
  geom_hline(yintercept = 0.5, linetype="dotdash",color = "grey")+
  xlab("log10(miRNA_Reads)")+
  ylab("miRNA_Genome_Ratio")+
  theme_bw()+
  theme(plot.margin = unit(c(1, 1, 1, 1),"cm"),
        legend.position=NaN,
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(face="bold", color="black", size=24,angle = 0,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))
}

#heatmap for miRNA
{
  Differential <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/miRNA/output/miRNA_GIvsNC_edger_exact_20220216.txt",header = TRUE, row.names = 1, sep = "\t")
  Differential_TPM <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/miRNA/matrix/miRNA_for_GIvsNC_20220216.txt",header = TRUE, row.names = 1, sep = "\t")
  #Differential_TPM <- as.data.frame(cpm(as.matrix(Differential_TPM)))
  
  Differential_TPM <- apply(as.matrix(Differential_TPM), 2, function(x) x / sum(as.numeric(x)) * 10^6)
  Differential_TPM <- as.data.frame(Differential_TPM)
  #samples <- des$samples
  #group <- des$group
  #mat_raw <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/miRNA/matrix/miRNA_for_GIvsNC_20220216.txt",check.names = FALSE,header = TRUE, row.names = 1, sep = "\t")
  #des <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/miRNA/group/des_miRNA_GIvsNC_20220216.csv", header = TRUE, check.names=FALSE, sep=',')
  #i=1
  #mat=as.data.frame(array(dim=c(length(rownames(mat_raw)),1)))
  #while (i<=length(samples)) {
  #  temp <- mat_raw[,which(colnames(mat_raw)==samples[i])]
  #  temp <- as.data.frame(temp)
  #  colnames(temp) <- samples[i]
  #  mat <- cbind(mat,temp)
  #  i=i+1
  #}
  
  #mat <- mat[,-1]
  #rownames(mat) <- rownames(mat_raw)
  #y <- DGEList(counts=mat, samples=samples, group=group)     #规定输入格式，必须包含这几个参数
  #keep <- filterByExpr(y,group = group)
  #y <- y[keep,]
  #y <- calcNormFactors(y, method="TMM")                      #归一化处理 TMM for RNA, CNV and Methylation
  

  Differential$Gene <- rownames(Differential)
  Differential_TPM$Gene <- rownames(Differential_TPM)
  
  forheatmap <- left_join(Differential,Differential_TPM,by=c("Gene"="Gene"))
  rownames(forheatmap) <- forheatmap$Gene
  forheatmap <- forheatmap[,-which(colnames(forheatmap)=="Gene")]

  gini_Expression_STAD <- gini(t(forheatmap[,grep("STAD",colnames(forheatmap))]))
  gini_Expression_STAD <- as.data.frame(gini_Expression_STAD)
  colnames(gini_Expression_STAD) <- c("gini_index_STAD")
  gini_Expression_CRC <- gini(t(forheatmap[,grep("CRC",colnames(forheatmap))]))
  gini_Expression_CRC <- as.data.frame(gini_Expression_CRC)
  colnames(gini_Expression_CRC) <- c("gini_index_CRC")
  gini_Expression_GI <- gini(t(forheatmap[,grep("STAD|CRC",colnames(forheatmap))]))
  gini_Expression_GI <- as.data.frame(gini_Expression_GI)
  colnames(gini_Expression_GI) <- c("gini_index_GI")
  gini_Expression_NC <- gini(t(forheatmap[,grep("NC",colnames(forheatmap))]))
  gini_Expression_NC <- as.data.frame(gini_Expression_NC)
  colnames(gini_Expression_NC) <- c("gini_index_NC")
  gini_Expression <- cbind(gini_Expression_CRC,gini_Expression_STAD,gini_Expression_GI,gini_Expression_NC)
  
  forheatmap <- cbind(gini_Expression, forheatmap)
  #forheatmap <- forheatmap[,-which(colnames(forheatmap)=="NC.PKU.mix0.qia")]
  #forheatmap <- forheatmap[,-which(colnames(forheatmap)=="NC.PKU.mix15.qia")]
  #forheatmap <- forheatmap[,-which(colnames(forheatmap)=="CRC.PKU.14.qia")]
  
  heatmap <- forheatmap[((forheatmap$pvalue < 0.05) & (forheatmap$log2FoldChange < 0)) | ((forheatmap$pvalue < 0.05) & (forheatmap$log2FoldChange > 0.67))
                        #& (forheatmap$gini_index_CRC < 0.5) & (forheatmap$gini_index_STAD < 0.5)
                        & (forheatmap$gini_index_GI < 0.5),]
                        # & (forheatmap$gini_index_NC < 1),]
  heatmap <- heatmap[order(heatmap$log2FoldChange),]
  rownames(heatmap) <- as.character(lapply(strsplit(rownames(heatmap),"\\."), function(x) x[1]))
  col_annotation <- data.frame(row.names = colnames(heatmap),CancerType=as.character(lapply(strsplit(colnames(heatmap),"\\."), function(x) x[1])))
  col_annotation$CancerType <- gsub("NC","HD",col_annotation$CancerType)
  col_annotation$CancerType <- factor(col_annotation$CancerType, levels = c("CRC","STAD","HD"), ordered = TRUE)
  heatmap_CRC <- heatmap[,grep("CRC",colnames(heatmap))]
  heatmap_STAD <- heatmap[,grep("STAD",colnames(heatmap))]
  heatmap_NC <- heatmap[,grep("NC",colnames(heatmap))]
  heatmap <- cbind(heatmap_CRC,heatmap_STAD,heatmap_NC)
  heatmap <- na.omit(heatmap)
  
  #library(biomaRt)
  #mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  
  value <- as.character(lapply(strsplit(rownames(heatmap),"\\."),function(x) x[1]))
  
  gene_names <- getBM(attributes=c("ensembl_transcript_id", "hgnc_symbol"),
                      filters = "ensembl_transcript_id",
                      values=value, mart= mart,useCache = FALSE)
  
  #gene_names[which(gene_names$ensembl_gene_id=="ENSG00000223911"),which(colnames(gene_names)=="hgnc_symbol")] <- "ENSG00000223911"
  #gene_names[which(gene_names$ensembl_gene_id=="ENSG00000235298"),which(colnames(gene_names)=="hgnc_symbol")] <- "ENSG00000235298"
  
  heatmap$ensembl_transcript_id <- as.character(lapply(strsplit(rownames(heatmap),"\\."),function(x) x[1]))
  heatmap <- left_join(heatmap,gene_names,by=c("ensembl_transcript_id"="ensembl_transcript_id"))
  rownames(heatmap) <- heatmap$hgnc_symbol
  heatmap <- heatmap[,-which(colnames(heatmap)=="ensembl_transcript_id")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="hgnc_symbol")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="gini_index_NC")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="gini_index_CRC")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="gini_index_STAD")]
  
  bk = unique(c(seq(-3,3, length=100)))
  pheatmap(
    heatmap,
    breaks = bk,
    annotation_col = col_annotation,
    #annotation_row = row_annotation,
    scale = "row",
    cluster_cols = FALSE,cluster_rows = FALSE,
    show_colnames= TRUE, 
    show_rownames=TRUE,
    #cluster_cols = FALSE,
    colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 10)
}

{
  Differential <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/miRNA/output/miRNA_CRCvsNC_edger_exact_20220216.txt",header = TRUE, row.names = 1, sep = "\t")
  Differential_TPM <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/miRNA/matrix/miRNA_for_CRCvsNC_20220216.txt",header = TRUE, row.names = 1, sep = "\t")
  #Differential_TPM <- as.data.frame(cpm(as.matrix(Differential_TPM)))
  
  Differential_TPM <- apply(as.matrix(Differential_TPM), 2, function(x) x / sum(as.numeric(x)) * 10^6)
  Differential_TPM <- as.data.frame(Differential_TPM)
  
  Differential$Gene <- rownames(Differential)
  Differential_TPM$Gene <- rownames(Differential_TPM)
  
  forheatmap <- left_join(Differential,Differential_TPM,by=c("Gene"="Gene"))
  rownames(forheatmap) <- forheatmap$Gene
  forheatmap <- forheatmap[,-which(colnames(forheatmap)=="Gene")]
  
  gini_Expression_CRC <- gini(t(forheatmap[,grep("CRC",colnames(forheatmap))]))
  gini_Expression_CRC <- as.data.frame(gini_Expression_CRC)
  colnames(gini_Expression_CRC) <- c("gini_index_CRC")
  gini_Expression_GI <- gini(t(forheatmap[,grep("STAD|CRC",colnames(forheatmap))]))
  gini_Expression_GI <- as.data.frame(gini_Expression_GI)
  colnames(gini_Expression_GI) <- c("gini_index_GI")
  gini_Expression_NC <- gini(t(forheatmap[,grep("NC",colnames(forheatmap))]))
  gini_Expression_NC <- as.data.frame(gini_Expression_NC)
  colnames(gini_Expression_NC) <- c("gini_index_NC")
  gini_Expression <- cbind(gini_Expression_CRC,gini_Expression_GI,gini_Expression_NC)
  
  forheatmap <- cbind(gini_Expression, forheatmap)
  #forheatmap <- forheatmap[,-which(colnames(forheatmap)=="NC.PKU.mix0.qia")]
  #forheatmap <- forheatmap[,-which(colnames(forheatmap)=="NC.PKU.mix15.qia")]
  #forheatmap <- forheatmap[,-which(colnames(forheatmap)=="CRC.PKU.14.qia")]
  
  heatmap <- forheatmap[((forheatmap$pvalue < 0.05) & (forheatmap$log2FoldChange < 0)) | ((forheatmap$pvalue < 0.05) & (forheatmap$log2FoldChange > 0))
                        #& (forheatmap$gini_index_CRC < 0.5) & (forheatmap$gini_index_STAD < 0.5)
                        & (forheatmap$gini_index_GI < 0.5),]
  # & (forheatmap$gini_index_NC < 1),]
  heatmap <- heatmap[order(heatmap$log2FoldChange),]
  rownames(heatmap) <- as.character(lapply(strsplit(rownames(heatmap),"\\."), function(x) x[1]))
  col_annotation <- data.frame(row.names = colnames(heatmap),CancerType=as.character(lapply(strsplit(colnames(heatmap),"\\."), function(x) x[1])))
  col_annotation$CancerType <- gsub("NC","HD",col_annotation$CancerType)
  col_annotation$CancerType <- factor(col_annotation$CancerType, levels = c("CRC","STAD","HD"), ordered = TRUE)
  heatmap_CRC <- heatmap[,grep("CRC",colnames(heatmap))]
  heatmap_NC <- heatmap[,grep("NC",colnames(heatmap))]
  heatmap <- cbind(heatmap_CRC,heatmap_NC)
  heatmap <- na.omit(heatmap)
  
  #library(biomaRt)
  #mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  
  value <- as.character(lapply(strsplit(rownames(heatmap),"\\."),function(x) x[1]))
  
  gene_names <- getBM(attributes=c("ensembl_transcript_id", "hgnc_symbol"),
                      filters = "ensembl_transcript_id",
                      values=value, mart= mart,useCache = FALSE)
  
  #gene_names[which(gene_names$ensembl_gene_id=="ENSG00000223911"),which(colnames(gene_names)=="hgnc_symbol")] <- "ENSG00000223911"
  #gene_names[which(gene_names$ensembl_gene_id=="ENSG00000235298"),which(colnames(gene_names)=="hgnc_symbol")] <- "ENSG00000235298"
  
  heatmap$ensembl_transcript_id <- as.character(lapply(strsplit(rownames(heatmap),"\\."),function(x) x[1]))
  heatmap <- left_join(heatmap,gene_names,by=c("ensembl_transcript_id"="ensembl_transcript_id"))
  rownames(heatmap) <- heatmap$hgnc_symbol
  heatmap <- heatmap[,-which(colnames(heatmap)=="ensembl_transcript_id")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="hgnc_symbol")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="gini_index_NC")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="gini_index_CRC")]
  
  bk = unique(c(seq(-3,3, length=100)))
  pheatmap(
    heatmap,
    breaks = bk,
    annotation_col = col_annotation,
    #annotation_row = row_annotation,
    scale = "row",
    cluster_cols = FALSE,cluster_rows = FALSE,
    show_colnames= TRUE, 
    show_rownames=TRUE,
    #cluster_cols = FALSE,
    colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 10)
}

{
  Differential <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/miRNA/output/miRNA_STADvsNC_edger_exact_20220216.txt",header = TRUE, row.names = 1, sep = "\t")
  Differential_TPM <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/miRNA/matrix/miRNA_for_STADvsNC_20220216.txt",header = TRUE, row.names = 1, sep = "\t")
  #Differential_TPM <- as.data.frame(cpm(as.matrix(Differential_TPM)))
  
  Differential_TPM <- apply(as.matrix(Differential_TPM), 2, function(x) x / sum(as.numeric(x)) * 10^6)
  Differential_TPM <- as.data.frame(Differential_TPM)
  
  Differential$Gene <- rownames(Differential)
  Differential_TPM$Gene <- rownames(Differential_TPM)
  
  forheatmap <- left_join(Differential,Differential_TPM,by=c("Gene"="Gene"))
  rownames(forheatmap) <- forheatmap$Gene
  forheatmap <- forheatmap[,-which(colnames(forheatmap)=="Gene")]
  
  gini_Expression_STAD <- gini(t(forheatmap[,grep("STAD",colnames(forheatmap))]))
  gini_Expression_STAD <- as.data.frame(gini_Expression_STAD)
  colnames(gini_Expression_STAD) <- c("gini_index_STAD")
  gini_Expression_GI <- gini(t(forheatmap[,grep("STAD|CRC",colnames(forheatmap))]))
  gini_Expression_GI <- as.data.frame(gini_Expression_GI)
  colnames(gini_Expression_GI) <- c("gini_index_GI")
  gini_Expression_NC <- gini(t(forheatmap[,grep("NC",colnames(forheatmap))]))
  gini_Expression_NC <- as.data.frame(gini_Expression_NC)
  colnames(gini_Expression_NC) <- c("gini_index_NC")
  gini_Expression <- cbind(gini_Expression_STAD,gini_Expression_GI,gini_Expression_NC)
  
  forheatmap <- cbind(gini_Expression, forheatmap)
  
  heatmap <- forheatmap[((forheatmap$padj < 0.05) & (forheatmap$log2FoldChange < 0.67)) | ((forheatmap$padj < 0.05) & (forheatmap$log2FoldChange > 0.67))
                        #& (forheatmap$gini_index_CRC < 0.5) & (forheatmap$gini_index_STAD < 0.5)
                        & (forheatmap$gini_index_GI < 0.5) & (forheatmap$gini_index_NC < 0.5),]
  heatmap <- heatmap[order(heatmap$log2FoldChange),]
  rownames(heatmap) <- as.character(lapply(strsplit(rownames(heatmap),"\\."), function(x) x[1]))
  col_annotation <- data.frame(row.names = colnames(heatmap),CancerType=as.character(lapply(strsplit(colnames(heatmap),"\\."), function(x) x[1])))
  col_annotation$CancerType <- gsub("NC","HD",col_annotation$CancerType)
  col_annotation$CancerType <- factor(col_annotation$CancerType, levels = c("CRC","STAD","HD"), ordered = TRUE)
  col_annotation <- na.omit(col_annotation)
  
  row_annotation <- data.frame(row.names = rownames(heatmap),log2FoldChange = heatmap$log2FoldChange)
  row_annotation$Trend <- ""
  row_annotation[row_annotation$log2FoldChange > 0,]$Trend <- "Up regulated in cancer"
  row_annotation[row_annotation$log2FoldChange < 0,]$Trend <- "Down regulated in cancer"
  row_annotation <- data.frame(row.names = rownames(row_annotation),Trend = row_annotation$Trend)
  row_annotation$Trend <- factor(row_annotation$Trend, levels = c("Down regulated in cancer","Up regulated in cancer"), ordered = TRUE)

  
  heatmap_STAD <- heatmap[,grep("STAD",colnames(heatmap))]
  heatmap_NC <- heatmap[,grep("NC",colnames(heatmap))]
  heatmap <- cbind(heatmap_STAD,heatmap_NC)
  heatmap <- na.omit(heatmap)
  
  #library(biomaRt)
  #mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  
  value <- as.character(lapply(strsplit(rownames(heatmap),"\\."),function(x) x[1]))
  
  gene_names <- getBM(attributes=c("ensembl_transcript_id", "hgnc_symbol"),
                      filters = "ensembl_transcript_id",
                      values=value, mart= mart,useCache = FALSE)
  
  
  heatmap$ensembl_transcript_id <- as.character(lapply(strsplit(rownames(heatmap),"\\."),function(x) x[1]))
  heatmap <- left_join(heatmap,gene_names,by=c("ensembl_transcript_id"="ensembl_transcript_id"))
  rownames(heatmap) <- heatmap$hgnc_symbol
  heatmap <- heatmap[,-which(colnames(heatmap)=="ensembl_transcript_id")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="hgnc_symbol")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="gini_index_NC")]
  heatmap <- heatmap[,-which(colnames(heatmap)=="gini_index_STAD")]
  
  row_annotation$ensembl_transcript_id <- rownames(row_annotation)
  row_annotation <- left_join(row_annotation,gene_names,by=c("ensembl_transcript_id"="ensembl_transcript_id"))
  row_annotation <- data.frame(row.names = row_annotation$hgnc_symbol,Trend = row_annotation$Trend)
  row_annotation$Trend <- factor(row_annotation$Trend, levels = c("Down regulated in cancer","Up regulated in cancer"), ordered = TRUE)
  
  bk = unique(c(seq(-3,3, length=100)))
  ann_color <- list(CancerType=c("CRC"="#EE7621","STAD"="red","HD"="blue"),Trend=c("Down regulated in cancer"="#D6D6D6","Up regulated in cancer"="#EB5E66"))
  ComplexHeatmap::pheatmap(
    heatmap,
    breaks = bk,
    annotation_col = col_annotation,
    annotation_row = row_annotation,
    annotation_colors=ann_color,
    gaps_row=c(26),
    gaps_col=c(13),
    scale = "row",
    cluster_cols = FALSE,cluster_rows = FALSE,
    show_colnames= FALSE, 
    show_rownames=TRUE,
    #cluster_cols = FALSE,
    colorRampPalette(c("black","#EEE9E9","red"))(100),
    fontsize_row = 10,
    name = "CPM")
}

#copy numebr plot
{
  library(ggplot2)
  library(ggpubr)
  bed_CRC <- read.table("/Users/yuhuan/Library/Containers/com.tencent.xinWeChat/Data/Library/Application Support/com.tencent.xinWeChat/2.0b4.0.9/25ce5949b711be7299e88400100eb330/Message/MessageTemp/81f849e734698b67ea205283f7c0b98e/File/未命名文件夹/CRC-PKU-10-wgs_bins.bed",sep = "\t",header = TRUE)
  bed_STAD <- read.table("/Users/yuhuan/Library/Containers/com.tencent.xinWeChat/Data/Library/Application Support/com.tencent.xinWeChat/2.0b4.0.9/25ce5949b711be7299e88400100eb330/Message/MessageTemp/81f849e734698b67ea205283f7c0b98e/File/未命名文件夹/STAD-PKU-11-wgs_bins.bed",sep = "\t",header = TRUE)
  bed_NC <- read.table("/Users/yuhuan/Library/Containers/com.tencent.xinWeChat/Data/Library/Application Support/com.tencent.xinWeChat/2.0b4.0.9/25ce5949b711be7299e88400100eb330/Message/MessageTemp/81f849e734698b67ea205283f7c0b98e/File/未命名文件夹/NC-PKU-10-wgs_bins.bed",sep = "\t",header = TRUE)
  plot_CRC <- bed_CRC[which(bed$chr=="1"),]
  plot_NC <- bed_NC[which(bed$chr=="1"),]
  plot_STAD <- bed_STAD[which(bed$chr=="1"),]
  
  p1 <- ggplot(plot_CRC,aes(x=start,y=ratio))+
    geom_point()+
    theme_bw()
  p2 <- ggplot(plot_STAD,aes(x=start,y=ratio))+
    geom_point()+
    theme_bw()
  p3 <- ggplot(plot_NC,aes(x=start,y=ratio))+
    geom_point()+
    theme_bw()
  
  ggarrange(p1,p2,p3,ncol = 1,align = "v")
}

#find genes
{
  library(biomaRt)
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  gene <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol","chromosome_name","start_position","end_position"),
                      filters = c("chromosome_name","start","end"),
                      values=list(chromosome="1",start="190000000",end="200000000"), mart= mart,useCache = FALSE)
}




