#Figure1 sample cohort
{
  ##multi-omics sample shaozhen
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/01.Sample summary/Figure3/")   
  
  #sunburst
  p <- read.csv("Figure3_forplot.csv",header = F, stringsAsFactors = FALSE)
  c <- read.csv("Color.csv",header = F,stringsAsFactors=FALSE)
  c[1,1] <- "NA"
  sunburst(p,
           count=TRUE,
           sortFunction = htmlwidgets::JS(
             "function(a,b) {
                       // sort by count descending
                       // unlike the other example using data.name, value is at the top level of the object
                       return b.value - a.value}" ),
           legend=list(w=240,h=30,r=10,s=5),
           #sortFunction = htmlwidgets::JS("Health","Colorectum","Liver","Stomach","Lung","Thyroid","Esophagus"),
           legendOrder = list("NC","CRC","STAD","Female","Male","Stage I","Stage II","Stage III","Stage IV","Unknown","long RNA","DNA methylation","WGS","small RNA"),
           colors=list(range = c$V2, domain = c$V1)
  )  
  
  #Age_Stage barplot
  clinical <- read.csv("./Figure3_subtypes_forplot.csv",header = TRUE)
  View(clinical)
  clinical[clinical==""]<-NA
  
  
  clinical$Stage.summary <- factor(clinical$Stage.summary,levels = c("missing","Stage I","Stage II","Stage III","Stage IV"))
  ggplot(data=clinical[which(clinical$Type!="NC"),],aes(x=Age_10,y=number,fill=Stage.summary))+ geom_bar(stat = "identity", width=0.9, col='transparent') +
    theme_bw()+
    theme(#legend.position="bottom",
      legend.position="right",
      panel.grid.major.x = element_blank(),
      legend.title = element_text(face="bold", color="transparent",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(color="black", size=20,angle = 0,hjust = 0.5),
      axis.text.y = element_text(color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    ylab("Number of patients")+xlab("Ages(years)")+
    scale_y_continuous(expand = c(0,0),limits = c(0,35))+
    scale_fill_manual(values=c("white","#DCDCDC","#A3A3A3","#4A4A4A","#000000"))
  
  
  
  #subtype
  clinical <- read.csv("./Figure3_simple_subtype_forplot.csv",header = TRUE,stringsAsFactors = FALSE)
  color <- read.csv("Color_for_subtype.csv",header = FALSE)
  clinical[clinical==""]<-NA
  clinical[is.na(clinical)] <- "No biopsy"
  
  color$V1 <- factor(color$V1,levels = color$V1)
  
  clinical$index <- factor(clinical$index,levels = c("CRC position","MMR","STAD position","HER2"))
  clinical$Subtype <- factor(clinical$Subtype,levels = color$V1)
  ggplot(data=clinical,aes(x=number,y=number,fill=clinical$Subtype))+ geom_bar(stat = "identity", width=10, col='transparent') + facet_grid(~index)+
    #geom_text(aes(label=Subtype), position = "stack",vjust=1)+
    theme_bw()+
    theme(
      legend.position="none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.title = element_text(face="bold", color="transparent",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(color="black", size=20,angle = 0,hjust = 0.5),
      axis.text.y = element_text(color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24),
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=10))+
    ylab("Number of patients")+xlab("")+scale_y_reverse()+scale_x_discrete(position = "top")+
    scale_fill_manual(values=as.character(color$V2))
  
}