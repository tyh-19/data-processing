```R
#gene coverage between RNA and DNA
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Expression/Coverage/")
DNA <- read.csv("GAPDH_DNA_NC.txt",sep = "\t",header = FALSE)
RNA <- read.csv("GAPDH_RNA_NC.txt",sep = "\t",header = FALSE)
colnames(DNA) <- c("Chr","Position","Depth")
colnames(RNA) <- c("Chr","Position","Depth")

RNA_DNA <- full_join(RNA,DNA,by = c("Position"="Position"))

colnames(RNA_DNA) <- c("Chr.RNA","Position","Depth.RNA","Chr.DNA","Depth.DNA") 
RNA_DNA <- RNA_DNA[,-4]
RNA_DNA <- RNA_DNA[,-1]

RNA_DNA[is.na(RNA_DNA)] <- 0

RNA_forplot <- RNA_DNA[,-3]
RNA_forplot$type <- "RNA"
colnames(RNA_forplot) <- c("Position","Depth","type")
DNA_forplot <- RNA_DNA[,-2]
DNA_forplot$type <- "DNA"
colnames(DNA_forplot) <- c("Position","Depth","type")
plot <- rbind(RNA_forplot,DNA_forplot)

plot$type <- factor(plot$type,levels=c("RNA","DNA"))
ggplot(plot,aes(x=plot$Position,y=log10(plot$Depth+1),fill=plot$type))+
  geom_bar(stat="identity",position = 'dodge')+
  scale_y_continuous(limits=c(0, log10(2000)),breaks = c(0,1,log10(2000)))+
  #scale_x_discrete(expand=c(0.2, 0.2))+
  #scale_fill_jco(alpha = 0.8)+
  scale_fill_manual(values=c("#5CACEE","#FFCC11"))+
  theme_bw()+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(size=1, colour = "black"),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=10,angle = 0,hjust = 1,vjust=0),
    axis.text.y = element_text(face="bold",  color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  xlab("")+
  ylab("Depth")+
  labs(fill="Type")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
gr <- GRanges("chr12", IRanges(6533927, 6538374), strand="+")
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                         org.Hs.eg.db,
                         gr=gr)

viewerStyle <- trackViewerStyle()
setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
trackList <- trackList(trs)
vp <- viewTracks(trackList, 
                 gr=gr, viewerStyle=viewerStyle, 
                 autoOptimizeStyle=TRUE)
addGuideLine(c(122929767, 122929969), vp=vp)
addArrowMark(list(x=122929650, 
                  y=1), # 2 means track 2 from the bottom.
             label="label",
             col="blue",
             vp=vp)
```

