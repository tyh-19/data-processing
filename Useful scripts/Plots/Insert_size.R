###insert length by miso
##length density
library(reshape2)
library(ggplot2)
plot <- read.csv("miso_Insert_length_summary.txt",sep = "\t",header = T)
head(plot)
plot <- plot[-1,]

colnames(plot)
#c(rep("Gastric Cancer",15),rep("Healthy Donor",15))
#colnames(plot) = c(rep("STAD",15),rep("Healthy_Donor",15))
#colnames(plot) = c(rep("sample",30))
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

{
  # freqency line plot
  summary(plot[,1]) #看数据的分布区间
  # set bin, min, max
  bin=1
  min=1
  max=1000
  m <- seq(min-1,max,by=bin)#设置一个区间范围
  barplot(table(cut(plot[,1],m)))#简单作图查看分布
  
  
  # 计算frequency
  freq = {}
  col <- colnames(plot) 
  i=1
  while(i<=length(col)){
    #prop.table(table(cut(plot$col[i],m)))*100#计算各个区间的频率
    freq.tmp <- as.data.frame(prop.table(table(cut(plot[,col[i]],m)))*100)
    freq <- cbind(freq,freq.tmp$Freq)
    i=i+1
  }
  # make frequency dataframe
  freq <- as.data.frame(freq)
  colnames(freq) <- colnames(plot)
  rownames(freq) <- as.data.frame(prop.table(table(cut(plot[,1],m))))[,1]
  
  # summary group frequency
  summary= {}
  summary$mean <- apply(freq, 1, mean,na.rm=T)
  #summary$mean_STAD <- apply(freq[,grep("STAD",colnames(freq))], 1, mean,na.rm=T)
  #summary$mean_HCC <- apply(freq[,grep("HCC",colnames(freq))], 1, mean,na.rm=T)
  #summary$mean_CRC <- apply(freq[,grep("CRC",colnames(freq))], 1, mean,na.rm=T)
  #summary$mean_ESCA <- apply(freq[,grep("ESCA",colnames(freq))], 1, mean,na.rm=T)
  #summary$mean_LUAD <- apply(freq[,grep("LUAD",colnames(freq))], 1, mean,na.rm=T)
  #summary$mean_NC <- apply(freq[,grep("NC",colnames(freq))], 1, mean,na.rm=T)
  summary.df <- as.data.frame(summary)
  summary.df$length <- seq(from=min-1,to=max-bin,by=bin)
  #write.csv(summary.df,"exoRBase_insert.csv")
  summary.df <- read.csv("plot_samtools_insert_3.csv",header=T)
  
  #plot
  ggplot(summary.df, aes(x=length,y=mean,colour=group)) + geom_line(size = 2) + theme_bw(base_size = 12, base_family = "Arial")+
    geom_vline(aes(xintercept = 167),linetype="dashed")+
    scale_x_continuous(limits=c(0, 1000),breaks = c(0,100,200,400,600,800,1000))+
    xlim(0,500)+
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
    labs(x="Insert Size (nt)",y="Freuqency (%)",title="", face="bold")
}



##insert length by RSEQC
plot <- read.csv("length_summary_Intronspanning.txt",sep = "\t",header = T,row.names=1)
plot <- plot[-201,]
freq = {}
col <- colnames(plot) 
i=1
while(i<=length(col)){
  #apply(as.data.frame(plot[,col[i]]), 2, sum, na.rm=T)
  freq.tmp <- as.data.frame(plot[,col[i]]/as.numeric(apply(as.data.frame(plot[,col[i]]), 2, sum, na.rm=T))*100)
  colnames(freq.tmp) <- "Freq"
  freq <- cbind(freq,freq.tmp$Freq)
  i=i+1
}

summary= {}

summary$mean <- apply(freq[,1:ncol(plot)], 1, mean,na.rm=T)
summary.df <- as.data.frame(summary)
summary.df$length <- summary.df$length <- seq(from=0,to=1000-5,by=5)

ggplot(summary.df, aes(x=length,y=mean_CRC)) + geom_line(color="dark blue",size = 2) + theme_bw(base_size = 12, base_family = "Arial")+
  geom_vline(aes(xintercept = 167),linetype="dashed")+
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
  labs(x="Insert Size (nt)",y="Freuqency (%)",title="", face="bold")


###insert length by samtools
##length density
library(reshape2)
library(ggplot2)

{
  ##if input file from Summary_insert_length.sh sep="\t"; if input file subset by awk '{NF=20}1', sep=" "
  plot <- read.csv("pico_length_summary_samtools.txt",sep = "\t",header = T)
  tail(plot)
  colnames(plot)
  # freqency line plot
  summary(plot[,1]) #看数据的分布区间
  # set bin, min, max
  bin=1
  min=1
  max=1000
  m <- seq(min-1,max,by=bin)#设置一个区间范围
  barplot(table(cut(plot[,1],m)))#简单作图查看分布
  
  
  # 计算frequency
  freq = {}
  col <- colnames(plot) 
  i=1
  while(i<=length(col)){
    #prop.table(table(cut(plot$col[i],m)))*100#计算各个区间的频率
    freq.tmp <- as.data.frame(prop.table(table(cut(plot[,col[i]],m)))*100)
    freq <- cbind(freq,freq.tmp$Freq)
    i=i+1
  }
  # make frequency dataframe
  freq <- as.data.frame(freq)
  colnames(freq) <- colnames(plot)
  rownames(freq) <- as.data.frame(prop.table(table(cut(plot[,1],m))))[,1]
  
  # summary group frequency
  summary= {}
  summary$mean <- apply(freq, 1, mean,na.rm=T)
  summary.df <- as.data.frame(summary)
  summary.df$length <- seq(from=min-1,to=max-bin,by=bin)
  write.csv(summary.df,"pico_freq_samtools.csv")
  
  #plot
  ggplot(summary.df, aes(x=length,y=mean)) + geom_line(color="dark blue",size = 2) + theme_bw(base_size = 12, base_family = "Arial")+
    geom_vline(aes(xintercept = 167),linetype="dashed")+
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
    labs(x="Insert Size (nt)",y="Freuqency (%)",title="", face="bold")
}