##PCA
library(ggplot2)
library(knitr)
library(psych)
library(openxlsx)
library(ggfortify)
pico_all_ID <- read.xlsx("C:/Users/Tao/Desktop/pico-all-ID.xlsx",rowNames=TRUE)
pico_all_for_pca <- read.csv("C:/Users/Tao/Desktop/pico-sample-for-pca.csv",header=TRUE)
pico_all_ID <- pico_all_ID[,1:4]
pico_all_ID
head(pico_all_for_pca)
##repeat without library_yield_log2
##princomp
pico_all_ID_PCA <- princomp(pico_all_ID,cor=TRUE)
summary(pico_all_ID_PCA,loadings=TRUE)
predict(pico_all_ID_PCA)
screeplot(pico_all_ID_PCA,type="lines")
biplot(pico_all_ID_PCA,choices=1:2,scale=1,pc.biplot=FALSE)
scoresdata=pico_all_ID_PCA$scores
head(scoresdata)
##prcomp
pico_all_ID_ggplot <- prcomp(pico_all_ID, scale = TRUE)
head(pico_all_ID_ggplot,1)
head(pico_all_ID_PCA)
##princomp and prcomp are the same
class(pico_all_ID_PCA)
pico_all_pcs <- data.frame(pico_all_ID_PCA$scores,Disease=pico_all_for_pca$Disease,Operator=pico_all_for_pca$Operator,Gender=pico_all_for_pca$Gender,Source=pico_all_for_pca$Source,Stage_summary=pico_all_for_pca$Stage_summary,Library_yield_log2=pico_all_for_pca$Library_yield_log2,Age=pico_all_for_pca$Age,Input_plasma=pico_all_for_pca$Input_plasma.ml.,RNA_mass=pico_all_for_pca$R._mass.ng.,PCR_cycles=pico_all_for_pca$PCR_cycles)
pico_all_pcs$PCR_cycles <- as.factor(pico_all_pcs$PCR_cycles)
head(pico_all_pcs)
write.csv(pico_all_pcs,"C:/Users/Tao/Desktop/pico-all_pcs_withoutlibraryyield.csv")

##percentage
percentage<-round(pico_all_ID_PCA$sdev/sum(pico_all_ID_PCA$sdev) * 100,2)
percentage<-paste(colnames(pico_all_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

##ggplot
ggplot(pico_all_pcs,aes(x=Comp.1,y=Comp.2))+
geom_point(aes(colour=Disease,shape=PCR_cycles),position="dodge")+
xlab(percentage[1])+ylab(percentage[2])+
theme_bw()+
#theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line= element_line(colour = "black"))+
facet_wrap(~Source)

#annotate('text', label = '14', x = -2, y = -1.25, size = 5, colour = '#f8766d')
#stat_ellipse(level = 0.95)+

##contribution of each feature
pico_all_ID_ggplot$rotation
PCA_r <- as.data.frame(pico_all_ID_ggplot$rotation)
PCA_r
PCA_r$feature <- row.names(PCA_r)
ggplot(PCA_r,aes(x=PC1,y=PC2,label=feature,color=feature )) + geom_point()+ geom_text(size=3,hjust=0.6,vjust=-0.7) + theme_classic()

##psych
fa.parallel(pico_all_ID, fa = 'pc', n.iter = 100, show.legend = FALSE)
pico_all_ID_psych <- principal(pico_all_ID, nfactors = 2, rotate = "none")
pico_all_ID_psych_score <- principal(pico_all_ID, nfactors = 2, rotate = "varimax", scores = TRUE)
options(digits = 2)
head(pico_all_ID_psych_score$scores)
autoplot(prcomp(pico_all_ID, scale = TRUE)) ##visualize by ggfortify

# 20200806 metagene
setwd("C:/Users/Tao/Desktop/Lulab/2020/02.Project/02. Data/metagene/DESeq2 and phyloseq/Unmapped_reads_classification")
pico_all_ID <- read.csv("counts-G.csv",header=TRUE)
rownames(pico_all_ID) <- pico_all_ID[,1]
pico_all_ID <- pico_all_ID[,-1]
head(pico_all_ID)
## to numeric in df by lapply
pico_all_ID=as.data.frame(lapply(pico_all_ID,as.numeric))
##princomp
pico_all_ID_PCA <- princomp(pico_all_ID,cor=TRUE)
summary(pico_all_ID_PCA,loadings=TRUE)
predict(pico_all_ID_PCA)
screeplot(pico_all_ID_PCA,type="lines")
biplot(pico_all_ID_PCA,choices=1:2,scale=1,pc.biplot=FALSE)
scoresdata=pico_all_ID_PCA$scores
head(scoresdata)
##prcomp
pico_all_ID
class(t(pico_all_ID))
pico_all_ID_ggplot <- prcomp(t(pico_all_ID), scale = TRUE)
head(pico_all_ID_ggplot)
length(pico_all_ID_ggplot)
class(pico_all_ID_ggplot)
##princomp and prcomp are the same
pico_all_for_pca <- read.csv("genus_for_plot.csv",header = TRUE)
class(pico_all_ID_PCA)
pico_all_pcs <- t(data.frame(pico_all_ID_PCA$scores,CANCERtypes=pico_all_for_pca$CANCERtypes,Health=pico_all_for_pca$Health))
head(pico_all_pcs)
write.csv(pico_all_pcs,"C:/Users/Tao/Desktop/pico-all_pcs_withoutlibraryyield.csv")

##percentage
percentage<-round(pico_all_ID_PCA$sdev/sum(pico_all_ID_PCA$sdev) * 100,2)
percentage<-paste(colnames(pico_all_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

##ggplot
ggplot(pico_all_pcs,aes(x=Comp.1,y=Comp.2))+
  geom_point(aes(colour=Disease,shape=PCR_cycles),position="dodge")+
  xlab(percentage[1])+ylab(percentage[2])+
  theme_bw()+
  #theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line= element_line(colour = "black"))+
  facet_wrap(~Source)

#annotate('text', label = '14', x = -2, y = -1.25, size = 5, colour = '#f8766d')
#stat_ellipse(level = 0.95)+

##contribution of each feature
pico_all_ID_ggplot$rotation
PCA_r <- as.data.frame(pico_all_ID_ggplot$rotation)
PCA_r
PCA_r$feature <- row.names(PCA_r)
ggplot(PCA_r,aes(x=PC1,y=PC2,label=feature,color=feature )) + geom_point()+ geom_text(size=3,hjust=0.6,vjust=-0.7) + theme_classic()

##psych
fa.parallel(pico_all_ID, fa = 'pc', n.iter = 100, show.legend = FALSE)
pico_all_ID_psych <- principal(pico_all_ID, nfactors = 2, rotate = "none")
pico_all_ID_psych_score <- principal(pico_all_ID, nfactors = 2, rotate = "varimax", scores = TRUE)
options(digits = 2)
head(pico_all_ID_psych_score$scores)
autoplot(prcomp(pico_all_ID, scale = TRUE)) ##visualize by ggfortify
