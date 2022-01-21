getwd()
##heatmap
library(pheatmap)
library(extrafont)
fonts()
TPM <- read.csv("TPM_filtered.txt",sep = '\t', header = T, row.names = 1)
head(TPM)
annotation_col = data.frame(CancerType=c("HD","HD","THCA","THCA","THCA","THCA","THCA","THN","THN"),row.names = c("Huaxi.NC.01","Huaxi.NC.02","THCA.0000687100","THCA.0002618932","THCA.0007244757","THCA.0008463257","THCA.0032983054","THN.0009164651","THN.0015322245"))
annotation_col
pheatmap(TPM, scale = "row", show_rownames = F, show_colnames = T, angle_col = 45, clustering_method = "complete", annotation_col = annotation_col)



##volcano
library(ggplot2)
library(ggrepel)
volcano <- read.csv("GC-NC-15v15.txt", sep = '\t', header = T, row.names = 1, fileEncoding = 'utf-8', fill = T)
head(volcano)
#volcano <- volcano[complete.cases(volcano), ]
#head(volcano)
#rownames(volcano) <- volcano[,1]
volcano$threshold <- as.factor(ifelse(volcano$padj < 0.05 & abs(volcano$log2FoldChange) >=1,ifelse(volcano$log2FoldChange > 1 ,'Up','Down'),'Not'))
write.csv(volcano,"THCA_volcano.csv")
smoke <- ggplot(data=volcano, aes(x=log2FoldChange, y =-log10(pvalue), color=threshold,fill=threshold)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_point(size=3) +
  #xlim(c(-6,8.5)) +
  #ylim(c(0,5)) +
  theme_bw(base_size = 12, base_family = "Arial") +
  geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6,aes(0,4))+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6,aes(-3,6))+
  theme(#legend.position="right",
    legend.position="none",
    panel.grid=element_blank(),
    #legend.title = element_blank(),
    #legend.text= element_text(face="bold", color="black",family = "Times", size=8),
    plot.title = element_text(hjust = 0.5,size=30,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=30),
    axis.text.y = element_text(face="bold",  color="black", size=30),
    axis.title.x = element_text(face="bold", color="black", size=36),
    axis.title.y = element_text(face="bold",color="black", size=36))+
  labs(x="log2(fold change)",y="-log10 (p value)",title="Different Gene Expression", face="bold")
smoke   

label <- read.csv("STAD_volcano.txt",sep="\t",header=T)
head(label,20)
smoke+geom_text_repel(aes(x=label$log2FoldChange, y =-log10(label$pvalue), label = label$X, vjust=0,hjust=0.4),size=8,color="black")
  #geom_text(aes(4,10,label="THCA"),hjust=0,vjust=0,colour="black",size=8)+geom_text(aes(-5,10,label="Healthy Donors"),hjust=0.1,vjust=0,colour="black",size=8)

##length density
library(reshape2)
library(ggplot2)
plot <- read.csv("length_summary.txt",sep = "\t",header = T)
head(plot)
plot <- plot[-1,]

colnames(plot)
c(rep("Gastric Cancer",15),rep("Healthy Donor",15))
colnames(plot) = c(rep("STAD",15),rep("Healthy_Donor",15))
colnames(plot) = c(rep("sample",30))
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

# freqency line plot
summary(plot$X.Huaxi.NC.01.hg38.remove.duplicates) #看数据的分布区间
# set bin, min, max
bin=10
min=1
max=400
m <- seq(min-1,max,by=bin)#设置一个区间范围
barplot(table(cut(plot$X.Huaxi.NC.01.hg38.remove.duplicates,m)))#简单作图查看分布


# 计算frequency
freq = {}
col <- colnames(plot) 
i=1
while(i<=length(col)){
  #prop.table(table(cut(plot$col[i],m)))*100#计算各个区间的频数
  freq.tmp <- as.data.frame(prop.table(table(cut(plot[,col[i]],m)))*100)
  freq <- cbind(freq,freq.tmp$Freq)
  i=i+1
}
# make frequency dataframe
freq <- as.data.frame(freq)
colnames(freq) <- colnames(plot)
rownames(freq) <- as.data.frame(prop.table(table(cut(plot$X.Huaxi.NC.01.hg38.remove.duplicates,m))))[,1]

# summary group frequency
summary= {}
summary$mean_STAD <- apply(freq[,1:15], 1, mean,na.rm=T)
summary$mean_NC <- apply(freq[,16:30], 1, mean,na.rm=T)
summary$mean <- apply(freq[,1:30], 1, mean,na.rm=T)
summary.df <- as.data.frame(summary)
summary.df$length <- seq(from=min,to=max-bin,by=bin)

#plot
ggplot(summary.df, aes(x=length,y=mean)) + geom_line(color="dark blue",size = 2) + theme_bw(base_size = 12, base_family = "Arial")+
  geom_vline(aes(xintercept = 166),linetype="dashed")+
  scale_x_continuous(limits=c(0, 1000),breaks = c(0,100,200,400,600,800,1000))+
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
  labs(x="Fragment Size (nt)",y="Freuqency (%)",title="", face="bold")




## RNA type



