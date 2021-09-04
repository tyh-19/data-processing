#prepare library and function
{
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(extrafont)
fonts()

qPCR_reshape <- function(data){
  # -△Ct =-(Candidate Ct - ref Ct)
  candidate <- rownames(qPCR[-c(1,2),])
  type <- unique(as.data.frame(t(qPCR[1,]))[,1])
  j=1
  qPCR.m={}
  while(j<=length(candidate)){
    target <- candidate[j]
    i=1
    gene={}
    while(i<=length(type)){
      temp <- qPCR[grep(target,rownames(qPCR),fixed=TRUE)[1],grep(as.character(type[i]),as.matrix(qPCR[1,]),fixed=TRUE)]
      hospital <- qPCR[2,grep(as.character(type[i]),as.matrix(qPCR[1,]),fixed=TRUE)]
      temp.t <- t(rbind(temp,hospital))
      colnames(temp.t) <- c("neg_deltaCt","Hospital")
      temp.df <- as.data.frame(temp.t)
      #k=1
      #OutVals <- boxplot(temp.t,plot=FALSE)$out
      #OutVals <- unique(OutVals)
      #OutVals <- sort(OutVals,decreasing = T)
      #while(k<=length(OutVals)){
      #  # return integer(0), temp.df <- temp.df
      #  if(length(grep(OutVals[k],temp.df$neg_deltaCt))==0){
      #    temp.df <- temp.df
      #  } else {
      #    temp.df <- as.data.frame(temp.df[-grep(OutVals[k],temp.df$neg_deltaCt),])
      #    colnames(temp.df) <- "neg_deltaCt"
      #    #print(length(temp.df$neg_deltaCt))
      #  }
      #  k=k+1
      #}
      temp.df$Group<-rep(type[i],times=length(temp.df$neg_deltaCt))
      temp.df$Candidates <- rep(candidate[j],times=length(temp.df$neg_deltaCt))
      gene <- rbind(gene,temp.df)
      i=i+1
    }
    qPCR.m <- rbind(qPCR.m,gene)
    j=j+1
  }
  
  qPCR.m$Candidates <- factor(qPCR.m$Candidates, levels = candidate)
  qPCR.m$neg_deltaCt <- as.numeric(as.character(qPCR.m$neg_deltaCt))
  qPCR.m <- unite(qPCR.m,"Candidates_Group",Candidates,Group,remove=FALSE)
  qPCR.m
}
}

#Remember to set work directory by setwd("Your/work/dir")

{
#load file
#file format: csv, colnames=sample_ID, 1st row=label by group, 2nd row=label by hosiptal, rownames=candidate_genes
#sample ID is unique for each sample
#group means two cohort for comparison, like exposure/non-exposure, cancer/nc
  ##sequncial in input file defines box sequence
  ##manually set group level to control box colour
{
qPCR <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/07.Pancancer_Biomarker/20210307/CRC_RPS8_RPL7A_RPL13_8vs15.csv", header = T, row.names = 1)
head(qPCR)
}

#reshape file
{
qPCR.m <- qPCR_reshape(qPCR)
print("These are input groups:")
print(as.character(unique(qPCR.m$Group)))
group <- factor(unique(qPCR.m$Group),levels=unique(qPCR.m$Group))
}

#plot
#base
{
qPCR.m$Candidates_Group <- factor(qPCR.m$Candidates_Group,levels=unique(qPCR.m$Candidates_Group))
qPCR.m$Group <- factor(qPCR.m$Group,levels=group)
p <- ggplot(qPCR.m,aes(x=Candidates_Group,y=neg_deltaCt,colour=Group, width = 5),fullrange = FALSE)+
  geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
  #geom_point(size = 3, position = position_jitterdodge(dodge.width=1,jitter.width = 0.1))+
  geom_point(aes(shape=Hospital),size = 3, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
  theme_bw()+
  theme(plot.margin = unit(c(1, 1, 3, 1),"cm"),
    legend.position="right",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(size=1, colour = "black"),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=30),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=30),
    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    axis.text.x = element_blank(),
    #axis.text.x = element_text(face="bold", color="black", size=30,angle = 45,hjust = 1),
    axis.text.y = element_text(face="bold",  color="black", size=30),
    axis.title.x = element_text(face="bold", color="black", size=36),
    axis.title.y = element_text(face="bold",color="black", size=36))+
  coord_cartesian(clip="off")+#key for annotate candidate 
  #https://stackoverflow.com/questions/64654407/annotate-ggplot2-across-both-axis-text-keeps-changing-position
  labs(x="",y="-delta Ct",title="", face="bold")+scale_color_aaas()
#scale_color_npg() 或 color_palette("npg")都是不错的配色方案
}

#add line
{
group_num=length(unique(qPCR.m$Group))
#https://stackoverflow.com/questions/44317502/update-a-ggplot-using-a-for-loop-r

i=1
while(i+1<=length(unique(qPCR.m$Candidates))){
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

#wilcox text by group, label FC and x.axis
##test between each 2 columns, positive:1st column,negative:2n columm. When not sure, check groups.
{
group_num=length(unique(qPCR.m$Group))  
  
i=1
j=1
groups <- unique(qPCR.m$Candidates_Group)
xlabel <- unique(qPCR.m$Candidates)
while(i<=length(unique(qPCR.m$Candidates_Group))){
  FC <- mean(subset(qPCR.m,Candidates_Group==groups[i])$neg_deltaCt,na.rm = TRUE)-mean(subset(qPCR.m,Candidates_Group==groups[i+1])$neg_deltaCt,na.rm = TRUE)
  #https://www.jianshu.com/p/aab34be5f983
  k=1
  l <- list()
  while(k<group_num){
    l_tmp <- list(c(groups[i],groups[i+k]))
    l <- c(l,l_tmp)
    k=k+1
  }
  if(FC >= 0 || FC=="NaN"){
  p <- p+stat_compare_means(comparison=l,
                            label = "p.signif",
                            size=5,
                            method="wilcox.test",
                            method.args = list(alternative = "greater"))+
    annotate("text",size=5,x=i+0.5,y=y_up-7,
             label=paste0("FC:",round(FC,2),"\nWilcox\nGreater"),family = "Arial", fontface = "bold")+
    #解释annotate: https://blog.csdn.net/g_r_c/article/details/19673625
    annotate("text",size=8,x=i+group_num/2-1+0.5,y=-Inf,
             label=xlabel[j],family = "Arial", fontface = "bold",angle=45,hjust=1,vjust=1)
    #geom_text(aes_(x=i+0.5,y=10,label=paste0("FC:",round(FC,3))))
  } else {
    p <- p+stat_compare_means(comparison=l,
                              label = "p.signif",
                              size=5,
                              method="wilcox.test",
                              method.args = list(alternative = "less"))+
      annotate("text",size=5,x=i+0.5,y=y_up-1,
               label=paste0("FC:",round(FC,2),"\nWilcox\nLess"),family = "Arial", fontface = "bold")+
      annotate("text",size=8,x=i+group_num/2-1+0.5,y=-Inf,
               label=xlabel[j],family = "Arial", fontface = "bold",angle=45,hjust=1,vjust=1)
  }
  i=i+group_num ## group_num controls gene name's position
  j=j+1
}
}

ggsave(p,filename = "CRC_RPS8_RPL7A_RPL13_8vs15.pdf",path = "./", width = 15, height = 10)
}

