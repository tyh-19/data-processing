library(ggpubr)
library(ggplot2)
library(tidyverse)
library(extrafont)
fonts()

TPM <- read.csv("TPM_all_pico.txt",sep="\t",header = T, row.names = 1)

head(grep("ENSG00000150991.14",rownames(TPM),value = T),10)

candidate <- c("ENSG00000075624.13","ENSG00000111640.14","ENSG00000166710.17","HMBS","ENSG00000142541.16","ENSG00000073578.16","ENSG00000150991.14","ENSG00000164924.17")
j=1
while(j<=length(candidate)){
target <- candidate[j]
type <- c("STAD","HCC","CRC","ESCA","LUAD","NC")
i=1
gene={}
while(i<=length(type)){
temp <- TPM[grep(target,rownames(TPM)),grep(type[i],colnames(TPM))]
temp.t <- t(temp)
colnames(temp.t) <- c("TPM")
temp.df <- as.data.frame(temp.t)
length(temp.t)
temp.df$CancerType<-rep(type[i],times=length(temp.t))
gene <- rbind(gene,temp.df)
i=i+1
}

#gene <- read.csv("test.csv")
#gene.m <- melt(gene,na.rm = TRUE)
my_comparisons <- list(c("NC", "HCC"), c("NC", "CRC"),c("NC","STAD"),c("NC","LUAD"), c("NC", "ESCA"))
gene$CancerType <- factor(gene$CancerType,levels=c("NC","HCC","CRC","STAD","LUAD","ESCA"))
p <- ggplot(gene,aes(x=CancerType,y=TPM))+
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
  labs(x="",y="TPM",title=target, face="bold")+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons)

ggsave(p,filename = paste0(target,".pdf"))
j=j+1
}
