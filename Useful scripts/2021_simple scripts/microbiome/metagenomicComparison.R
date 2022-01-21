#### environment preparation
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(tidyverse)
library(extrafont)
library(edgeR)
library(DESeq2)
fonts()
#this is a package for outlier removal, seems hard to install
#devtools::install_github('kongdd/Ipaper')
#library(Ipaper)

#### differential analysis
### DESeq2
{
mat <- read.table("counts-G_filtered.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
des <- read.csv("des_CRCvsLUAD.csv", header = TRUE, check.names=FALSE, sep=',')
samples <- des$samples
group <- des$group
batch <- des$batch
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = as.matrix(des),
                              design = ~group+batch)
dds <- DESeq(dds)
res <- results(dds, contrast=c('group', 'positive', 'negative'))
write.table(as.data.frame(res), "deseq2_LUAD.txt", sep='\t', quote=FALSE, row.names=TRUE)
}

### EdgeR_glmlrt
{
mat <- read.table("counts-G_filtered.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
des <- read.csv("des_CRCvsLUAD.csv", header = TRUE, check.names=FALSE, sep=',')
samples <- des$samples
group <- des$group
batch <- des$batch
y <- DGEList(counts=mat, samples=samples, group=group)
y <- calcNormFactors(y, method="TMM")
design <- model.matrix(~group+batch)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
test <- glmLRT(fit, coef=2)
res <- topTags(test, n=nrow(mat), sort.by='none')
res <- cbind(res$table, baseMean=2^(res$table$logCPM))
mapped_names <- colnames(res)
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
write.table(res, "edger_glmlrt_LUAD.txt", sep='\t', quote=FALSE, row.names=TRUE)
}

#### correclation between deseq2 and edger_glmlrt
library(ggplot2)
{
deseq <- read.csv("deseq2_STAD.txt", sep = "\t", header =TRUE, row.names = 1)
egder <- read.csv("edger_glmlrt_STAD.txt", sep = "\t", header =TRUE, row.names = 1)

head(egder)

ggplot(deseq,aes(deseq$log2FoldChange)) +
  geom_histogram()

FC={}
FC$egder <- egder$log2FoldChange
FC$deseq <- deseq$log2FoldChange
FC <- as.data.frame(FC)
rownames(FC) <- row.names(deseq)

r <- cor(FC$egder,FC$deseq,method="pearson")
ggplot(data=FC, aes(x=egder, y=deseq)) +
  geom_point(color="#d7191c") +
  #geom_smooth(method="lm",color="#1a9641") +
  geom_text(aes(x=2, y=-3,label=paste("R","=",signif(r,3),seq="")),color="#fdae61",size=12)+
  theme_bw()+
  xlab("Edger")+
  ylab("DESeq2")+
  theme(
    axis.title = element_text(face="bold", color="black", size=24),
    axis.text = element_text(face="bold",  color="black", size=24)
  )
}
#### figure1 Differential expression between CRC and others
volcano <- read.csv("egder_glmlrt_decontamed.txt", sep = '\t', header = T, row.names = 1, fileEncoding = 'utf-8', fill = T)
volcano$threshold <- as.factor(ifelse(volcano$padj < 0.05 & abs(volcano$log2FoldChange) >=1,ifelse(volcano$log2FoldChange > 1 ,'Up','Down'),'Not'))
write.csv(volcano,"CRCvsOthers_volcano.csv")
abundant <- c("1350|Enterococcus","1883|Streptomyces", "1301|Streptococcus","32008|Burkholderia","2093|Mycoplasma",
              "1485|Clostridium","1279|Staphylococcus","2745|Halomonas","150247|Anoxybacillus","1912216|Cutibacterium",
              "196118|Methanocaldococcus","1716|Corynebacterium","662|Vibrio","745|Pasteurella","107|Spirosoma",
              "482|Neisseria","613|Serratia","1578|Lactobacillus","673534|Plantactinospora","561|Escherichia")
k=1
while(k<=length(abundant)){
volcano <- volcano[-grep(abundant[k],rownames(volcano)),]
k=k+1
}
head(volcano)
#volcano <- volcano[complete.cases(volcano), ]
#head(volcano)
#rownames(volcano) <- volcano[,1]
volcano_plot <- ggplot(data=volcano, aes(x=log2FoldChange, y =-log10(pvalue), colour=threshold)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_point(size=volcano$baseMean/sum(volcano$baseMean)*300) +
  #xlim(c(-6,8.5)) +
  #ylim(c(0,5)) +
  theme_bw(base_size = 12, base_family = "Arial") +
  geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6,aes(0,4))+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6,aes(-3,6))+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    #legend.title = element_blank(),
    #legend.text= element_text(face="bold", color="black",family = "Times", size=8),
    plot.title = element_text(hjust = 0.5,size=30,face="bold"),
    axis.text.x = element_text(face="bold",family = "Arial",color="black", size=30),
    axis.text.y = element_text(face="bold",family = "Arial",color="black", size=30),
    axis.title.x = element_text(face="bold",family = "Arial",color="black", size=36),
    axis.title.y = element_text(face="bold",family = "Arial",color="black", size=36))+
  labs(x="log2(fold change)",y="-log10 (p.adj)",title="Different Microbe Abundance", face="bold")
volcano_plot   



# need to modify CRCvsOthers_volcano.csv manually
label <- read.csv("CRCvsOthers_volcano.csv",header=T,row.names = 1)
k=1
while(k<=length(abundant)){
  label<- label[-grep(abundant[k],rownames(label)),]
  k=k+1
}
head(label,20)
volcano_plot+geom_text_repel(aes(x=label$log2FoldChange, y =-log10(label$pvalue), label = label$X, vjust=0,hjust=0.5),size=5,color="black")
#geom_text(aes(4,10,label="THCA"),hjust=0,vjust=0,colour="black",size=8)+geom_text(aes(-5,10,label="Healthy Donors"),hjust=0.1,vjust=0,colour="black",size=8)



#### boxplot to illustrate batch effect (raw abundance normalized to relative abundance)
## read in files
#raw to illustrate batch
counts <- read.csv("counts-G_filtered_norm.txt",sep="\t",header = T, row.names = 1)

#check distribution
head(grep("Cyprinivirus",rownames(counts),value = T),10)
distribution <- as.data.frame(t(counts))
ggplot(distribution,aes(distribution$`1279|Staphylococcus`)) +
  geom_histogram()

## input candidates
candidate <- c("Varicellovirus","Calothrix","Roseomonas","Thiohalobacter","Geobacillus",
               "Oleiphilus","Salegentibacter","Nodularia","Phreatobacter","Pseudothermotoga","Crenobacter",
               "Ictalurivirus","Marinithermus","Chamaesiphon","Kordia","Archaeoglobus","Belliella","Candidatus Portiera",
               "Raphidiopsis","Gillisia","Trichormus","Chondrocystis","Odoribacter","Azotobacter","Rivularia","Natranaerobius",
               "Halotalea","Roseolovirus","Aquitalea","Candidatus Nitrosocaldus","Rhodoluna","Hyperthermus","Chroococcidiopsis",
               "Vitreoscilla","Geminocystis","Plantactinospora","Marichromatium","Jeongeupia")
genome_size <- c(1,6.96039,6.56385,3.29962,4.6871,6.56089,4.21216,5.77354,4.97505,
                 2.14743,3.42805,1,2.26917,1.47113,5.48239,2.70174,4.64319,0.411975,
                 3.18651,4.37796,7.01825,5.00252,5.87234,5.08446,8.72877,3.19145,4.38743,
                 0.173861,4.47281,1.58767,1.42451,0.969948,6.91578,2.61267,
                 4.42606,9.21918,3.76327,3.78881)

j=1

### figure1. compare all microbe on one boxplot
{
  j=1
  total={}
  while(j<=length(candidate)){
    target <- candidate[j]
    type <- c("CRC","HCC","STAD","LUAD","ESCA","NC")
    i=1
    gene={}
    while(i<=length(type)){
      temp <- counts[grep(target,rownames(counts)),grep(type[i],colnames(counts))]
      temp.t <- t(temp)
      temp.t <- temp.t/genome_size[j]
      colnames(temp.t) <- c("counts")
      temp.df <- as.data.frame(temp.t)
      k=1
      OutVals = boxplot(temp.t,plot=FALSE)$out
      while(k<=length(OutVals)){
        temp.df <- as.data.frame(temp.df[-grep(OutVals[k],temp.df$counts),])
        colnames(temp.df) <- "counts"
        k=k+1
      }
      temp.df$CancerType<-rep(type[i],times=length(temp.df$counts))
      temp.df$microbe<-rep(target,times=length(temp.df$counts))
      gene <- rbind(gene,temp.df)
      i=i+1
    }
    total <- rbind(total,gene)
    j=j+1
  }
  
  total$CancerType <- factor(total$CancerType,levels=type)
  total$microbe <- factor(total$microbe,level=candidate)
  
  ggplot(total,aes(x=microbe,y=counts,fill=CancerType))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1),outlier.size=0,outlier.alpha = 0)+
    #geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 0.5))+
    #ylim(0,0.3)+
    scale_fill_manual(values=c("#EE0000", "#F0F8FF","#B9D3EE", "#4372AA", "#26466D", "#308014")) + 
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
      axis.title.y = element_text(face="bold",color="black", size=24))+
    #geom_boxplot()+
    stat_compare_means(aes(group = CancerType),label = "p.signif")+
    labs(x="",y="Relative Abundance (%) normalized by genome size",title="CRC microbes", face="bold")
}

## figure1.1 compare all microbe on one boxplot with breaks
down <- ggplot(total,aes(x=microbe,y=counts,fill=CancerType))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1),outlier.size=0,outlier.alpha = 0)+
  #geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
  scale_fill_manual(values=c("#EE0000", "#F0F8FF", "#B9D3EE", "#4372AA", "#26466D", "#308014")) +
  theme_bw()+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(size=1, colour = "black"),
    legend.title = element_blank(),
    legend.text= element_blank(),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
    axis.text.y = element_text(face="bold",  color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  labs(x="",y="",title="", face="bold")+
  #geom_boxplot()+
  #stat_compare_means(comparisons = my_comparisons)+
  coord_cartesian(ylim = c(0,0.01))+
  scale_y_continuous(breaks = c(0, 0.005,0.01))

middle <- ggplot(total,aes(x=microbe,y=counts,fill=CancerType))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1),outlier.size=0,outlier.alpha = 0)+
  #geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
  scale_fill_manual(values=c("#EE0000", "#F0F8FF", "#B9D3EE", "#4372AA", "#26466D", "#308014")) +
  theme_bw()+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(size=1, colour = "black"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(face="bold",  color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  #geom_boxplot()+
  #stat_compare_means(comparisons = my_comparisons)+
  labs(x="",y="Relative Abundance(%) normalized by genome size",title="", face="bold")+
  coord_cartesian(ylim = c(0.02, 0.1))+
  scale_y_continuous(breaks = c(0.02,0.06,0.1))

upper <- ggplot(total,aes(x=microbe,y=counts,fill=CancerType))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1),outlier.size=0,outlier.alpha = 0)+
  #geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
  scale_fill_manual(values=c("#EE0000", "#F0F8FF", "#B9D3EE", "#4372AA", "#26466D", "#308014")) +
  theme_bw()+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(size=1, colour = "black"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(face="bold",  color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  #geom_boxplot()+
  #stat_compare_means(comparisons = my_comparisons)+
  labs(x="",y="",title="", face="bold")+
  coord_cartesian(ylim = c(0.7, 10))+
  scale_y_continuous(breaks = c(0.7, 6, 9))

ggarrange(upper,middle,down,ncol = 1,nrow = 3,align = "v",common.legend = TRUE)



#### figure2.1 boxplot for each microbe to illustrate batch effect
j=1
{
while(j<=length(candidate)){
target <- candidate[j]
type <- c("STAD","HCC","CRC","ESCA","LUAD","NC")
i=1
gene={}
while(i<=length(type)){
temp <- counts[grep(target,rownames(counts)),grep(type[i],colnames(counts))]
temp.t <- t(temp)
colnames(temp.t) <- c("counts")
temp.df <- as.data.frame(temp.t)
k=1
OutVals = boxplot(temp.t,plot=FALSE)$out
while(k<=length(OutVals)){
  temp.df <- as.data.frame(temp.df[-grep(OutVals[k],temp.df$counts),])
  colnames(temp.df) <- "counts"
  k=k+1
}
temp.df$CancerType<-rep(type[i],times=length(temp.df$counts))
temp.df$microbe<-rep(target,times=length(temp.df$counts))
gene <- rbind(gene,temp.df)
i=i+1
}

#gene <- read.csv("test.csv")
#gene.m <- melt(gene,na.rm = TRUE)
my_comparisons <- list(c("CRC","STAD"),c("CRC", "HCC"),c("CRC","LUAD"), c("CRC", "ESCA"), c("CRC", "NC"))
gene$CancerType <- factor(gene$CancerType,levels=c("CRC","STAD","HCC","LUAD","ESCA","NC"))
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
    axis.title.y = element_text(face="bold",color="black", size=24))+
  #geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons)+
  labs(x="",y="Relative Abundance (%)",title=target, face="bold")
  
ggsave(p,filename = paste0(target,".pdf"))
j=j+1
}
}

## figure2.2 detailed boxplot for batch effect
# ESCA batch: ESCA_KZ,ESCA_PKU,ESCA_TJHH,ESCA_UnionH,ESCA_XieH
# HCC batch: HCC.ShHW,HCC_ChQ,HCC_ShH
# NC batch: NC_ChQ,NC_HaiB,NC_PKU
#raw to illustrate batch
counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/CRCvsOthers_considerbatch/matrix/counts-G_filtered_norm.txt",sep="\t",header = T, row.names = 1)
{
j=1
test={}
while(j<=length(all_microbes)){
  target <- all_microbes[j]
  type <- c("CRC","NC_PKU","NC_ShH","NC_ChQ","NC","HCC","HCC.ShHW","HCC_ChQ","HCC_ShH","STAD","LUAD","ESCA")
  i=1
  gene={}
  while(i<=length(type)){
    temp <- counts[grep(target,rownames(counts),fixed=TRUE),grep(type[i],colnames(counts),fixed=TRUE)]
    temp.t <- t(temp)
    colnames(temp.t) <- c("counts")
    temp.df <- as.data.frame(temp.t)
    #k=1
    #OutVals = boxplot(temp.t,plot=FALSE)$out
    #while(k<=length(OutVals)){
    #  temp.df <- as.data.frame(temp.df[-grep(OutVals[k],temp.df$counts),])
    #  colnames(temp.df) <- "counts"
    #  k=k+1
    #}
    temp.df$CancerType<-rep(type[i],times=length(temp.df$counts))
    gene <- rbind(gene,temp.df)
    i=i+1
  }

  #gene <- read.csv("test.csv")
  #gene.m <- melt(gene,na.rm = TRUE)
  my_comparisons <- list(c("CRC","NC_PKU"),c("CRC","NC_ShH"),c("CRC","NC_ChQ"),c("CRC","NC"),c("CRC","HCC"),c("CRC","HCC.ShHW"),c("CRC","HCC_ChQ"), c("CRC","HCC_ShH"), c("CRC", "STAD"),c("CRC","LUAD"),c("CRC","ESCA"))
  gene$CancerType <- factor(gene$CancerType,levels=type)
  p <- ggplot(gene,aes(x=CancerType,y=counts,fill=CancerType))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
    geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 0.6))+
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
      axis.title.y = element_text(face="bold",color="black", size=24))+
    #geom_boxplot()+
    stat_compare_means(comparisons = my_comparisons,
                       method = "wilcox.test",
                       method.args = list(alternative = "two.sided")
                       #label = "p.signif"
                       )+
    labs(x="",y="Relative Abundance (%)",title=target, face="bold")
  
  
  ggsave(p,filename = paste0(target,"_NC_batch.pdf"))
  t <- compare_means(counts ~ CancerType, gene, 
                     method = "wilcox.test", 
                     #not support one tailed test
                     #method.args = list(alternative = "greater"), 
                     p.adjust.method = "BH")
  t$microbe <- rep(target,times=length(t$group1))
  test <- rbind(test,t)
  j=j+1
}
  write.csv(test,"all_microbes_wilcoxon.two.sided.csv")
}

