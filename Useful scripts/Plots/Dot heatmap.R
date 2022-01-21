library(ggplot2)
library(plyr)
library(reshape2)
library(scales)
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2020/02.Project/01.Sequencing and Cancer/01. Cancer plasma microbiome/Plasma_Tissue/GSEA/20201028")
#import data
kreuz <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
kreuz.m <- melt(kreuz)
head(kreuz)
write.csv(kreuz,"kreuz.csv")
head(kreuz.m)
write.csv(kreuz.m,"kreuz.m.csv")

#reshape data
plot <- read.csv("FDR_1028.csv", header = T)
plot.m <- melt(plot, value.name = "FDR")
head(plot)
head(plot.m,10)
write.csv(plot.m,"FDR_1028.m.csv")

#import data
plot <- read.csv("ES_FDR_1028.m.revised.csv", header = T)
plot$GeneSet <- factor(plot$GeneSet,levels=c("COAD_Tumor_vs_Normal_Up_500","READ_Tumor_vs_Normal_Up_500","CRC_Tumor_vs_Normal_Up_500","COAD_Tumor_vs_Normal_Down_500","READ_Tumor_vs_Normal_Down_500","CRC_Tumor_vs_Normal_Down_500"))
plot$Rank <- factor(plot$Rank, levels=c("CRC.2418658.plasma_vs_average_HD_PKU","CRC.2418503.plasma_vs_average_HD_PKU","CRC.2418488.plasma_vs_average_HD_PKU","CRC.2416785.plasma_vs_average_HD_PKU","CRC.2418277.plasma_vs_average_HD_PKU","Tissue_average_Tumor_vs_average_Normal"))
p <- ggplot(plot, aes(GeneSet, Rank)) +
    geom_point(aes(size = FDR, colour = NES))+
    geom_text(aes(x=GeneSet,y=Rank,label=FDR,hjust=0.5,vjust=3), size=3,color="black",position = position_dodge(width=0.00),check_overlap = FALSE)+
    #geom_tile(aes(fill = value), colour = "white") +
    scale_size(breaks = waiver(), trans = "reverse",range = c(1, 10))+
    #scale_size_binned(breaks = waiver(),n.breaks = 10, trans = "reverse",range = c(1, 10))+
    scale_color_gradient2(breaks=waiver(), name="NES", low ="blue", mid= ("white"), high = "red", midpoint = 0)
base_size <- 10
p + theme_grey(base_size = base_size) +
  theme_bw()+
  labs(x = "GeneSet", y = "Rank", title = "GSEA")+
  theme(legend.position = "right", axis.ticks = element_blank(), 
        axis.text.x = element_text(size = base_size *0.8, angle = 45, hjust = 1,
                                   colour = "black", family = "Arial"),
        axis.text.y = element_text(size = base_size *0.8,
                                   colour = "black", family = "Arial"))

