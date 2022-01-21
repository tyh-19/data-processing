setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/")
library(openxlsx)
library(ggplot2)
library(extrafont)
fonts()
## sunburst
{
devtools::install_github("timelyportfolio/sunburstR")

library(sunburstR)
library(htmlwidgets)
# read in sample visit-sequences.csv data provided in source
# https://gist.github.com/kerryrodden/7090426#file-visit-sequences-csv
sequences <- read.csv(
  system.file("examples/visit-sequences.csv",package="sunburstR")
  ,header=F
  ,stringsAsFactors = FALSE
)

View(sequences)
sunburst(sequences)


p <- read.csv("RNA.csv",header = F, stringsAsFactors = FALSE)
sunburst(p,
         count=TRUE, 
         legend=list(w=240,h=30,r=10,s=5)
         #colors=list("#BFBFBF","#D3D3D3","#E0E0E0"),
         )

#how to set color: https://stackoverflow.com/questions/49993198/how-to-specify-the-colors-and-toggle-labels-for-each-category-in-r-sunburst
p <- read.csv("cfRNA.csv",header = F, stringsAsFactors = FALSE)

labels <- c("Health","Colorectum","Liver","Stomach","Lung","Thyroid","Esophagus","long RNA","small RNA","DNA methylation","DNA seq")
colors <- c("#8FD8D8","#FF3030","#FF7F50","#FCB514","#FFFFAA","#FEE5AC","#E8F1D4","#B0B0B0","#E0E0E0","#F0F0F0","#FCFCFC")


sunburst(p,
         #height = 10,
         #width = 20,
         count=TRUE, 
         legendOrder = list("Health","Colorectum","Liver","Stomach","Lung","Thyroid","Esophagus","long RNA","small RNA","DNA methylation","DNA seq"),
         legend=list(w=150,h=30,r=10,s=5)
         #colors=list(range = colors, domain = labels)
)

## Future plan
p <- read.csv("FuturePlan.csv",header = F, stringsAsFactors = FALSE)
labels <- c("Health","Colorectum","Liver","Stomach","Lung","Thyroid","Esophagus","RNA","DNA","DNA methylation","Protein+Metabolite","Early Stage","Other Stage","No Stage")
colors <- c("#8FD8D8","#FF3030","#FF7F50","#FCB514","#FFFFAA","#FEE5AC","#E8F1D4","#A3A3A3","#A8A8A8","#BDBDBD","#D6D6D6","#7AA9DD","#162252","#FFFFFF")


sunburst(p,
         #height = 10,
         #width = 20,
         count=TRUE, 
         sortFunction = htmlwidgets::JS(
           "function(a,b) {
            // sort by count descending
            // unlike the other example using data.name, value is at the top level of the object
            return b.value - a.value}" ),
         #sortFunction = htmlwidgets::JS("Health","Colorectum","Liver","Stomach","Lung","Thyroid","Esophagus"),
         #legendOrder = list("Health","Colorectum","Liver","Stomach","Lung","Thyroid","Esophagus","long RNA","small RNA","DNA methylation","DNA seq"),
         legend=list(w=150,h=30,r=10,s=5)
         #colors=list(range = colors, domain = labels)
)

##Published RNA data
{
p <- read.csv("Published_RNA.csv",header = F, stringsAsFactors = FALSE)
c <- read.csv("Color_published_RNA.csv",header = F)
sunburst(p,
         count=TRUE,
         sortFunction = htmlwidgets::JS(
                      "function(a,b) {
                       // sort by count descending
                       // unlike the other example using data.name, value is at the top level of the object
                       return b.value - a.value}" ),
         legend=list(w=240,h=30,r=10,s=5),
         colors=list(range = c$V2, domain = c$V1)
           )
}
}

##sample sumary
colors <- c("#8FD8D8","#FF3030","#FF7F50","#FCB514","#FFFFAA","#FEE5AC","#E8F1D4","#A3A3A3","#A8A8A8")
sunburst(p,
         #height = 10,
         #width = 20,
         count=TRUE, 
         sortFunction = htmlwidgets::JS(
           "function(a,b) {
            // sort by count descending
            // unlike the other example using data.name, value is at the top level of the object
            return b.value - a.value}" ),
         #sortFunction = htmlwidgets::JS("Health","Colorectum","Liver","Stomach","Lung","Thyroid","Esophagus"),
         #legendOrder = list("Health","Colorectum","Liver","Stomach","Lung","Thyroid","Esophagus","long RNA","small RNA","DNA methylation","DNA seq"),
         legend=list(w=150,h=30,r=10,s=5))
         #colors=list(range = colors, domain = labels)

##multi-omics sample shaozhen
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA")      
p <- read.csv("Figure1_forplot.csv",header = F, stringsAsFactors = FALSE)
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
         
         
         
data <- read.xlsx("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/exRNA-datasets-processed_status_revised.xlsx",sheet=1)

##barplot
{
nsum <- function(num){
  sum(-num)
} 

data$Disease_Type = with(data, reorder(Disease_Type, Sample_Size, nsum))
ggplot(data[-length(data$Disease_Type),],aes(x=Disease_Type,y=Sample_Size,fill=RNA_Source))+
  geom_col()+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(size=1, colour = "black"),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=18,angle = 45,hjust = 1),
    axis.text.y = element_text(face="bold",  color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  labs(x="Disease Type",y="Sample Size",title="RNA published data", face="bold")
}

##piechart
{
specimen <- unique(data[-length(data$Specimen),]$Specimen)
i=1
plot={}
while(i<=length(specimen)){
  tmp <- data[which(data$Specimen==specimen[i],arr.ind=TRUE),]
  pie <- data.frame(Sample_Size=0,Specimen=0)
  pie$Sample_Size <- apply(as.matrix(tmp$Sample_Size),2,sum)
  pie$Specimen <- specimen[i]
  plot <- rbind(plot,pie)
  i=i+1
}

ggplot(plot,aes(x="",y=Sample_Size,fill=Specimen))+
  geom_bar(stat="identity",alpha=1.6)+
  scale_fill_brewer(palette="Set1")+
  coord_polar(theta="y")+
  xlab("")+
  ylab("")+
  labs(title="RNA published data",fill="Specimen")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.grid=element_blank(),
        legend.position="right",
        panel.border=element_blank(),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        #axis.text.x = element_text(face="bold", color="black", size=18,angle = 45,hjust = 1),
        #axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
  geom_text(aes(x=1.6,y=10),label=paste("Total Sample Size: ",3428),size=8)
  geom_text(aes(y=cumsum((plot$Sample_Size)-plot$Sample_Size/2),x=1),label=plot$Sample_Size)
}
