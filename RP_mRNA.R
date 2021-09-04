#### environment preparation
{
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(tidyverse)
  library(reshape)
  library(extrafont)
  library(edgeR)
  library(DESeq2)
  fonts()
  library(progress)
  library(pheatmap)
  library(clusterProfiler) ## enrichment
  library(biomaRt) ## ID convert
  library(org.Hs.eg.db) ## annotation for human
  library(enrichplot) ##gseaplot
  library(GO.db)
  library(dplyr)
}
#functions
{
  ##TPM calculator (referred to ensembl ID longest transcript)
  TPM <- function(rawcount,output,transcript_length) {
    counts <- rawcount
    values <- as.character(lapply(strsplit(rownames(counts),".",fixed = TRUE),function(x) x[1]))
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    annotations <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"), 
                         filters = "ensembl_gene_id", 
                         values=unique(values), mart= mart,useCache = FALSE)
    library(dplyr)
    longest_transcript <- arrange(annotations, ensembl_gene_id, dplyr::desc(transcript_length))
    write.csv(longest_transcript,transcript_length)
    #merge 会改变顺序，尽量不要用
    #merged <- merge(as.data.frame(rownames(counts)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by.x=1,by.y="hgnc_symbol",all.x=TRUE)
    
    #left_join 不会改变顺序
    merged <- left_join(as.data.frame(values),longest_transcript[!duplicated(longest_transcript$ensembl_gene_id),],by = c("values"="ensembl_gene_id"))
    
    #对于无注释长度的转录本，记为1000，即不对长度标准化
    merged$transcript_length[is.na(merged$transcript_length)] <- 1000
    
    transcript_lengths <- as.vector(merged$transcript_length)
    
    #find gene length normalized values 
    rpk <- apply( counts, 2, function(x) x/(transcript_lengths/1000))
    #normalize by the sample size using rpk values
    tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
    
    new_rownames <- unite(merged, "new_rownames", values, hgnc_symbol, transcript_length, sep = "|")
    
    rownames(tpm) <- new_rownames$new_rownames
    tpm <- as.data.frame(tpm)
    rownames(tpm) <- gsub(".","|",fixed = TRUE,rownames(tpm))
    tpm$gene_id <- new_rownames$new_rownames
    tpm <- aggregate(. ~ gene_id, data = tpm, sum)
    rownames(tpm) <- tpm$gene_id
    tpm <- tpm[,-which(colnames(tpm)=="gene_id")]
    write.table(tpm,output,quote = FALSE,sep="\t")
  }
}
##workdir
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA")

#rawcounts to TPM
{
path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_rawcount/"
raw_matrix <- dir(path)
i=1
while(i<=length(raw_matrix)) {
prefix <- as.character(lapply(strsplit(raw_matrix[i],".",fixed = TRUE),function(x) x[1]))
rawcount <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_rawcount/",raw_matrix[i]),sep = "\t", header = TRUE, row.names = 1,check.names = FALSE)
output <- paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/",prefix,"_TPM.txt")
transcript_length <- paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/transcript_length/",prefix,"_longest_transcript.csv")
TPM(rawcount,output,transcript_length)
i=i+1
}
}

#get single pathway gene TPM mean
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA")
path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/"
TPM_matrix <- dir(path)

prefix <- as.character(lapply(strsplit(TPM_matrix,"_TPM",fixed = TRUE),function(x) x[1]))

i=1
while(i<=length(TPM_matrix)){
counts <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/",TPM_matrix[i]),sep = "\t", header = TRUE,check.names = FALSE)
#pathway_gene <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/gene/RP-RNA.csv",header = T)
pathway_gene <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/gene/Ribosome.csv",header = T)

#get pathway genes
j=1
pathway_gene_count={}
while(j<=nrow(pathway_gene)){
  target <- pathway_gene[j,2]
  gene_symbol <- pathway_gene[j,1]
  if(length(grep(target,rownames(counts)))==0) {
    #print(paste0("No ",target," in this dataset."))
    temp <- as.data.frame(array(,dim=c(1,ncol(counts))))
    temp[1,] <- 0
    rownames(temp) <- gene_symbol
    j=j+1
  } else {
    #temp <- counts[which(rownames(counts)==target),]
    temp <- counts[grep(target,rownames(counts),fixed=TRUE),]  #for ensg
    rownames(temp) <- gene_symbol
    pathway_gene_count <- rbind(pathway_gene_count,temp)
    j=j+1
  }
}
pathway_gene_count_log2 <- log2(pathway_gene_count+1)
pathway_gene_count_log2_colmean <- colMeans(pathway_gene_count_log2)
pathway_gene_count_log2_colmean <- as.data.frame(pathway_gene_count_log2_colmean)
colnames(pathway_gene_count_log2_colmean) <- "Ribosome"
pathway_gene_count_colmean <- colMeans(pathway_gene_count)
pathway_gene_count_colmean <- as.data.frame(pathway_gene_count_colmean)
colnames(pathway_gene_count_colmean) <- "Ribosome"

dir.create(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/Ribosome_matrix/",prefix[i]))
write.table(pathway_gene_count,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/Ribosome_matrix/",prefix[i],"/",prefix[i],"_Ribosome_matrix.txt"),quote = FALSE,sep = "\t")
write.table(pathway_gene_count_log2_colmean,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/Ribosome_matrix/",prefix[i],"/",prefix[i],"_Ribosome_log2.txt"),quote = FALSE,sep = "\t")
write.table(pathway_gene_count_colmean,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/Ribosome_matrix/",prefix[i],"/",prefix[i],"_Ribosome.txt"),quote = FALSE,sep = "\t")
i=i+1
}
}

#RP-mRNA ratio in Ribosome
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA")
path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/"
TPM_matrix <- dir(path)

prefix <- as.character(lapply(strsplit(TPM_matrix,"_TPM",fixed = TRUE),function(x) x[1]))

i=1
RP_mRNA_ratio <- {}
while(i<=length(prefix)){
Ribosome <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/Ribosome_matrix/",prefix[i],"/",prefix[i],"_Ribosome.txt"),header = TRUE, row.names = 1,sep = "\t")
RP_mRNA <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/RP_mRNA_matrix/",prefix[i],"/",prefix[i],"_RP-RNA.txt"),header = TRUE, row.names = 1,sep = "\t")

RP_mRNA_ratio_tmp <- mean((RP_mRNA$RP_mRNA*78)/(Ribosome$Ribosome*182))
RP_mRNA_ratio <- list(RP_mRNA_ratio,RP_mRNA_ratio_tmp)
i=i+1
}
RP_mRNA_ratio <- as.data.frame(unlist(RP_mRNA_ratio))
RP_mRNA_ratio$dataset <- prefix
colnames(RP_mRNA_ratio) <- c("RP_mRNA_ratio","dataset")

ggplot(RP_mRNA_ratio[-c(2,3,7,9,11,18),],aes(x=dataset,y=RP_mRNA_ratio,fill=RP_mRNA_ratio))+
  geom_bar(stat = "identity",colour = "black")+
  #geom_text(aes(x=name,y=performance-0.85+0.005,label=performance),size = 4,angle = 0)+
  #scale_y_continuous(breaks = c(0,0.05,0.10,0.15),labels = c("0.85","0.90","0.95","1"),expand = c(0,0),limits = c(0,0.15))+
  #scale_fill_manual(values = RNA_color)+
  xlab("")+
  ylab("RP_mRNA ratio")+
  theme_bw()+
  theme(#legend.position="right",
    plot.margin = unit(x=c(15,5,10,5),units="pt"),
    legend.position="null",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_text(face="bold", color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
    axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
    axis.title.x = element_text(face="bold", color="black", size=16),
    axis.title.y = element_text(face="bold",color="black", size=20))
}

#RP-mRNA boxplot in tissue
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA")
path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/RP_mRNA_matrix/"
Pathway_mean <- grep("TCGA",dir(path),value = TRUE)

prefix <- as.character(lapply(strsplit(Pathway_mean,"_count_matrix",fixed = TRUE),function(x) x[1]))

i=1
while(i<=length(Pathway_mean)){
des <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/group_info/",prefix[i],"_PrimaryTumor_vs_TissueNormal.csv"),header = TRUE, row.names = 1)
log2_mean <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/RP_mRNA_matrix/",Pathway_mean[i],"/",Pathway_mean[i],"_RP_mRNA_log2.txt"),sep = "\t",header = TRUE, row.names = 1,check.names=FALSE)
rownames(log2_mean) <- gsub(".","-",fixed= TRUE,rownames(log2_mean))
cancer_tissue <- rownames(des[des$group=="positive",])
normal_tissue <- rownames(des[des$group=="negative",])
#des <- des[,-which(colnames(des)=="batch")]
des$sample_id <- rownames(des)
log2_mean$sample_id <- rownames(log2_mean)
boxplot <- left_join(log2_mean,des,by = c('sample_id'='sample_id'))

boxplot <- boxplot[,-which(colnames(boxplot)=="batch")]
boxplot <- boxplot[which(boxplot$group!="Recurrent Solid Tumor"),]
boxplot <- na.omit(boxplot)

my_comparisons <- list(c("positive","negative"))
p <- ggplot(boxplot,aes(x=group,y=RP_mRNA,fill=group))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
  geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
  scale_fill_manual(values = c("blue","red")) +
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

m=mean(boxplot[boxplot$group=="positive",]$RP_mRNA)-mean(boxplot[boxplot$group=="negative",]$RP_mRNA)
if(m>0){
  p <- p+stat_compare_means(comparisons = my_comparisons,
                            method = "wilcox.test",
                            method.args = list(alternative = "greater"),
                            label = "p.signif"
  )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.greater"), face="bold",fill="Type")
  ggsave(p,device = "pdf",path="/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/boxplot/",filename = paste0(prefix[i],"_greater.pdf"),width = 6,height = 8)
} else if(m==0){
  p <- p+stat_compare_means(comparisons = my_comparisons,
                            method = "wilcox.test",
                            method.args = list(alternative = "two.sided"),
                            label = "p.signif"
  )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.twosided"), face="bold",fill="Type")
  ggsave(p,device = "pdf",path="/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/boxplot/",filename = paste0(prefix[i],"_twosided.pdf"),width = 6,height = 8)
} else {
  p <- p+stat_compare_means(comparisons = my_comparisons,
                            method = "wilcox.test",
                            method.args = list(alternative = "less"),
                            label = "p.signif"
  )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.less"), face="bold",fill="Type")
  ggsave(p,device = "pdf",path="/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/boxplot/",filename = paste0(prefix[i],"_less.pdf"),width = 6,height = 8)
}
i=i+1
}
}

#RP-mRNA boxplot in blood
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA")
  path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/RP_mRNA_matrix/"
  Pathway_mean <- grep("TCGA",dir(path),value = TRUE,invert = TRUE)
  Pathway_mean <- grep("TEP",Pathway_mean,value = TRUE,invert = TRUE)
  
  prefix <- as.character(lapply(strsplit(Pathway_mean,"_",fixed = TRUE),function(x) x[1]))
  
  i=1
  while(i<=length(Pathway_mean)){
    des <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/group_info/",prefix[i],"_all_for_boxplot.csv"),header = TRUE, row.names = 1)
    log2_mean <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/RP_mRNA_matrix/",Pathway_mean[i],"/",Pathway_mean[i],"_RP_mRNA_log2.txt"),sep = "\t",header = TRUE, row.names = 1,check.names=FALSE)
    des$sample_id <- rownames(des)
    log2_mean$sample_id <- rownames(log2_mean)
    boxplot <- left_join(log2_mean,des,by = c('sample_id'='sample_id'))
    
    boxplot <- boxplot[,which(colnames(boxplot)!="batch")]
    boxplot <- boxplot[which(boxplot$group!="Recurrent Solid Tumor"),]
    boxplot <- boxplot[which(boxplot$group!="Other_disease"),]
    boxplot <- na.omit(boxplot)
    
    my_comparisons <- list(c("positive","negative"))
    p <- ggplot(boxplot,aes(x=group,y=RP_mRNA,fill=group))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
      scale_fill_manual(values = c("blue","red")) +
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
    
    m=mean(boxplot[boxplot$group=="positive",]$RP_mRNA)-mean(boxplot[boxplot$group=="negative",]$RP_mRNA)
    if(m>0){
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                method.args = list(alternative = "greater"),
                                label = "p.signif"
      )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.greater"), face="bold",fill="Type")
      ggsave(p,device = "pdf",path="/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/boxplot/",filename = paste0(prefix[i],"_greater.pdf"),width = 6,height = 8)
    } else if(m==0){
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                method.args = list(alternative = "two.sided"),
                                label = "p.signif"
      )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.twosided"), face="bold",fill="Type")
      ggsave(p,device = "pdf",path="/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/boxplot/",filename = paste0(prefix[i],"_twosided.pdf"),width = 6,height = 8)
    } else {
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                method.args = list(alternative = "less"),
                                label = "p.signif"
      )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.less"), face="bold",fill="Type")
      ggsave(p,device = "pdf",path="/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/boxplot/",filename = paste0(prefix[i],"_less.pdf"),width = 6,height = 8)
    }
    i=i+1
  }
}

#detailed boxplot
des <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/group_info/GSE89843_all.csv"),header = TRUE, row.names = 1)
log2_mean <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/RP_mRNA_matrix/GSE89843_featurecounts/GSE89843_featurecounts_RP_mRNA_log2.txt"),sep = "\t",header = TRUE, row.names = 1,check.names=FALSE)
des$sample_id <- rownames(des)
log2_mean$sample_id <- rownames(log2_mean)
boxplot <- left_join(log2_mean,des,by = c('sample_id'='sample_id'))

boxplot <- boxplot[,which(colnames(boxplot)!="batch")]
boxplot <- boxplot[which(boxplot$group!="Recurrent Solid Tumor"),]
boxplot <- boxplot[which(boxplot$group!="Other_disease"),]
boxplot <- boxplot[which(boxplot$group!="nonCancer"),]
boxplot <- na.omit(boxplot)

my_comparisons <- list(c("HD","NSCLC"),c("HD","Chronic Pancreatitis"),c("HD","Multiple Sclerosis"),
                       c("HD","Epilepsy"),c("HD","Pulmonary Hypertension"),c("HD","Non-significant Atherosclerosis"),
                       c("HD","Stable Angina Pectoris"),c("HD","Unstable Angina Pectoris"))
boxplot$group <- factor(boxplot$group,levels=c("HD","NSCLC","Chronic Pancreatitis","Multiple Sclerosis",
                                               "Epilepsy","Pulmonary Hypertension","Non-significant Atherosclerosis","Stable Angina Pectoris","Unstable Angina Pectoris"))
p <- ggplot(boxplot,aes(x=group,y=RP_mRNA,fill=group))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
  geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
  scale_fill_manual(values = c("blue","red","#EE7621","#EE7621","#EE7621","#EE7621","#EE7621","#EE7621","#EE7621")) +
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

#m=mean(boxplot[boxplot$group=="positive",]$RP_mRNA)-mean(boxplot[boxplot$group=="negative",]$RP_mRNA)
m=1
if(m>0){
  p <- p+stat_compare_means(comparisons = my_comparisons,
                            method = "wilcox.test",
                            method.args = list(alternative = "greater"),
                            label = "p.signif"
  )+labs(x="",y="log2(TPM+1)",title=paste0("GSE89843","\nwilcox.test.greater"), face="bold",fill="Type")
  #ggsave(p,device = "pdf",path="/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/boxplot/",filename = paste0(prefix[i],"_greater.pdf"),width = 6,height = 8)
} else if(m==0){
  p <- p+stat_compare_means(comparisons = my_comparisons,
                            method = "wilcox.test",
                            method.args = list(alternative = "two.sided"),
                            label = "p.signif"
  )+labs(x="",y="log2(TPM+1)",title=paste0("GSE89843","\nwilcox.test.twosided"), face="bold",fill="Type")
  #ggsave(p,device = "pdf",path="/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/boxplot/",filename = paste0(prefix[i],"_twosided.pdf"),width = 6,height = 8)
} else {
  p <- p+stat_compare_means(comparisons = my_comparisons,
                            method = "wilcox.test",
                            method.args = list(alternative = "less"),
                            label = "p.signif"
  )+labs(x="",y="log2(TPM+1)",title=paste0("GSE89843","\nwilcox.test.less"), face="bold",fill="Type")
  #ggsave(p,device = "pdf",path="/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/boxplot/",filename = paste0(prefix[i],"_less.pdf"),width = 6,height = 8)
}
