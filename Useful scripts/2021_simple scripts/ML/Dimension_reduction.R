setwd("/Users/yuhuan/Desktop/")
library(Rtsne)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggsci)

target <- "."

# data
counts <- read.csv("counts-G_norm.txt",sep="\t",header = T,row.names = 1)
type <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/04.plasma microbiome subtype and abundant microbes/dimension reduction/type.csv",header= T)

count.t <- as.data.frame(t(counts[,grep(target,colnames(counts))]))

## tSNE
{
set.seed(42)
##nrow(count.t) larger than 3*perplexity
tsne_out <- Rtsne(count.t,pca=FALSE,dims=2,perplexity = 30,theta=0.0) 

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)

## 拟杆菌Bacteroides、普氏菌Prevotella和瘤胃球菌Ruminococcus 
ggplot(tsne_res,aes(tSNE1,tSNE2,color=type[grep(target,type$taxo),]$Ruminococcus_norm)) + 
  geom_point() + theme_bw() + 
  scale_color_gradient(low = "blue",high = "red",) +
  #geom_hline(yintercept = 0,lty=2,col="red") + 
  #geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  stat_ellipse(aes(tSNE1,tSNE2, color = type[grep(target,type$taxo),]$Ruminococcus_norm), geom = "polygon", alpha = 0.3, levles = 0.99)+
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Ruminococcus",color="Type")
}


## PCA
PCA_out <- PCA(count.t, graph = FALSE)
fviz_eig(PCA_out, addlabels = TRUE, ylim = c(0, 50))
#eig.val结果中第一列是特征值，第二列是可解释变异的比例，第三列是累计可解释变异的比例
eig.val <- get_eigenvalue(PCA_out)


fviz_pca_ind(PCA_out,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = type[grep(target,type$taxo),]$type, # color by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Type"
)
