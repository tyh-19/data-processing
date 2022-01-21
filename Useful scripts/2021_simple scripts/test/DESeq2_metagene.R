library(DESeq2)
setwd("C:/Users/Tao/Desktop/DESeq2 and phyloseq/LUAD")
LUAD_countdata <- read.csv("CRC_Counts-S.csv", header = T)
LUAD_coldata <- read.table("CRC_group.txt",header = T)
head(LUAD_countdata)
head(LUAD_coldata)

rownames(LUAD_countdata) <- LUAD_countdata[,1]
LUAD_countdata <- LUAD_countdata[,-1]

LUAD_dds <- DESeqDataSetFromMatrix(LUAD_countdata, LUAD_coldata, design= ~ group)

LUAD_dds_normalized <- DESeq(LUAD_dds)

LUAD_res_normalized = results(LUAD_dds_normalized)

LUAD_res_normalized = LUAD_res_normalized[order(LUAD_res_normalized$pvalue),]

summary(LUAD_res_normalized)

head(LUAD_res_normalized)

table(LUAD_res_normalized$padj<0.05)

write.csv(LUAD_res_normalized,file="CRC_diff.csv")

##volcano
library(ggplot2)
library(ggrepel)
volcano <- read.csv("LUAD_smoke_diff_microbiome-COFG.csv",header=T)
volcano <- volcano[complete.cases(volcano), ]
head(volcano)
rownames(volcano) <- volcano[,1]
volcano$threshold <- as.factor(ifelse(volcano$padj < 0.05 & abs(volcano$log2FoldChange) >=1,ifelse(volcano$log2FoldChange > 1 ,'Up','Down'),'Not'))
write.csv(volcano,"LUAD_smoke_label-COFG_all.csv")
smoke <- ggplot(data=volcano, aes(x=log2FoldChange, y =-log10(padj), color=threshold,fill=threshold)) +
scale_color_manual(values=c("blue", "grey","red"))+
geom_point(size=3) +
xlim(c(-6,8.5)) +
ylim(c(0,5)) +
theme_bw(base_size = 12, base_family = "Times") +
geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6,aes(0,4))+
geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6,aes(-3,6))+
theme(#legend.position="right",
      legend.position="none",
      panel.grid=element_blank(),
      #legend.title = element_blank(),
      #legend.text= element_text(face="bold", color="black",family = "Times", size=8),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
      labs(x="log2(fold_change)",y="-log10 (q_value)",title="Different microbial genus of smoke condition in LUAD", face="bold")
smoke   

label <- read.csv("LUAD_smoke_label-COFG.csv",header=T)
head(label)
smoke+geom_text_repel(aes(x=label$log2FoldChange, y =-log10(label$padj), label = label$X, vjust=0.7,hjust=-0.1),size=5)+
   geom_text(aes(4,5,label="Never_Smokers"),hjust=0,vjust=0,colour="black",size=8)+geom_text(aes(-5,5,label="Current_Smokers"),hjust=0.1,vjust=0,colour="black",size=8)

#label <- subset(volcano,threshold!="Not",select=c("log2FoldChange","padj"))
#head(label)
#label$genus <- rownames(label)

#ggsave("/home/test/share/RNAseq_homework/volcano_plot_tyh_genelabe_repel.pdf", plot=volcanoplot_gene, height = 10, width = 10)
