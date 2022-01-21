#limma
#接着按照文档的说明以及limma包的习惯，我们需要对count进行标准化以及转化为log2的值，这里标准化的方法为TMM，使用edgeR里面的calcNormFactors函数即可
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)

GC_saliva_RPM <- read.table("C:/Users/Tao/Desktop/Candidates/STAD/相关数据/saliva GC biomarker/GSE121870_DWONG1-gastric-cancer-and-controls-2016-10-25_exceRpt_gencode_ReadsPerMillion.txt",header = TRUE)
GC_saliva_group<- read.csv("C:/Users/Tao/Desktop/Candidates/STAD/相关数据/saliva GC biomarker/GC_saliva_design.csv", header =TRUE)
head(GC_saliva_RPM)
head(GC_saliva_group)
log2GC_saliva_RPM <- log2(GC_saliva_RPM)
head(log2GC_saliva_RPM[1:10,1:3])

group_list <- factor(GC_saliva_group$Group)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(GC_saliva_RPM)

class(GC_saliva_RPM)
str(GC_saliva_RPM)
#limma-trend
fit <- lmFit(log2GC_saliva_RPM, design)
fit <- eBayes(fit, trend=TRUE)
plotSA(fit, main="Final model: Mean-variance trend")
output <- topTable(fit, coef=ncol(design), n=Inf)
output
sum(output$adj.P.Val<0.05)
str(output)
write.csv(output,"C:/Users/Tao/Desktop/Candidates/STAD/相关数据/saliva GC biomarker/GC_saliva_diffgenes.csv")


#volcano plot
#remove NA then plot
diffgene <- na.omit(output)
colnames(diffgene)
diffgene$threshold <- as.factor(ifelse(diffgene$adj.P.Val < 0.05 & abs(diffgene$logFC) >=1,ifelse(diffgene$logFC > 1 ,'Up','Down'),'Not'))
summary(diffgene$threshold)
ggplot(data=diffgene, aes(x=logFC, y =-log10(adj.P.Val), color=threshold,fill=threshold)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_point(size=1) +
  xlim(c(-7, 7)) +
  ylim(c(0,80)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        #legend.key.size = unit(4,'cm'),
        legend.text= element_text(face="bold", color="black", size=10),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=12),
        axis.text.y = element_text(face="bold",  color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(x="log2(fold_change)",y="-log10 (q_value)",title="Volcano plot of different expression gene", face="bold")

#
diffgene_FC2 <- diffgene[-log(diffgene$adj.P.Val, base=10) > 120,]
diffgene_FC2 <- diffgene_FC2[abs(diffgene$logFC) >= 1,]
str(diffgene_FC2)
rownames(diffgene_FC2)

library(DESeq2)
library(pheatmap)
diffgene_FC2_RPMmat <- GC_saliva_RPM[which(rownames(GC_saliva_RPM) %in% rownames(diffgene_FC2)),]
str(diffgene_FC2_RPMmat)
z_score_matrix <- t(scale(t(diffgene_FC2_RPMmat)))
annotation_col <- data.frame(GC_saliva_group[,3])
head(annotation_col)
rownames(annotation_col) <- GC_saliva_group[,1]
colnames(annotation_col) <- c("sample_type")
annotation_col

colnames(z_score_matrix) <- rownames(annotation_col)

pheatmap(mat = z_score_matrix, cluster_rows = T, show_colnames = T, cluster_cols = T,
         annotation_col = annotation_col,clustering_distance_cols = "correlation",
         annotation_colors = list(Sample_type = c(Cancer="orange", Control='blue')),
         show_rownames = T, color = colorRampPalette(c("navy","white","firebrick3"))(100),
         legend = T, legend_breaks = c(-2,0,2), breaks = unique(seq(-2, 2, length = 100)),
         cutree_cols = 2, cutree_rows = 2,
         legend_labels = c('-2','0','2'))

##venn
install.packages("VennDiagram")
library(VennDiagram)
inhouse_down <- read.table("C:/Users/Tao/Desktop/Candidates/STAD/相关数据/saliva GC biomarker/cfRNA_inhouse_down_1511.txt",header = F)
saliva_down <- read.table("C:/Users/Tao/Desktop/Candidates/STAD/相关数据/saliva GC biomarker/GC_saliva_down_330.txt",header = F)
inhouse_down <- inhouse_down$V1
saliva_down <- saliva_down$V1

length(inhouse_down); length(saliva_down)

venn.diagram(list(inhouse_down = inhouse_down, saliva_down = saliva_down), fill = c("yellow", "cyan"), cex = 1.5, filename = "C:/Users/Tao/Desktop/Candidates/STAD/相关数据/saliva GC biomarker/Sample1_VS_Sample2_down.png")

inter <- get.venn.partitions(list(inhouse_down = inhouse_down, saliva_down = saliva_down))

for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(3, 4)], 'C:/Users/Tao/Desktop/Candidates/STAD/相关数据/saliva GC biomarker/venn_inter_down.txt', row.names = FALSE, sep = '\t', quote = FALSE)

##GO
inter_up <- read.table("C:/Users/Tao/Desktop/Candidates/STAD/相关数据/saliva GC biomarker/up/up.txt",header = T,sep="\t",quote = "")
inter_up[10,]
colnames(inter_up)
ggplot(data=inter_up[1:10,])+
  geom_bar(aes(x=reorder(Term,Fold.Enrichment),y=Fold.Enrichment, fill=-log10(PValue)),
           stat='identity') +
  coord_flip() +
  scale_fill_gradient(expression(-log["10"](P.value)),low="blue", high = "red")+
  theme_light()+
  theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.text.y=element_text(family="myFont",face="bold",size=12),
    axis.ticks.length = unit(0, "cm"),
    #axis.ticks=element_blank(),
    #axis.line = element_line(colour="black",size=1.2,linetype="solid"),
    #legend.position = "bottom",
    legend.text = element_text(family="myFont",face="bold",size=14),
    legend.direction = NULL,
    legend.title = element_text(family="myFont",face="bold",size=14),
    #legend.key = element_rect(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill="white",color="transparent"),
    panel.border = element_rect(linetype="solid",,size=0.7,colour = "black"))+
  theme(title=element_text(family="myFont",face="bold",size=14),
        text=element_text(family="myFont",face="bold",size=38),
        axis.text=element_text(family="myFont",face="bold",size=10,colour = "black"))+
  xlab("") +
  ylab("Fold Enrichment")

str(plasma_genetype)
##piechart
library(ggplot2)
saliva_genetype <- read.csv("C:/Users/Tao/Desktop/Candidates/STAD/相关数据/saliva GC biomarker/saliva_gene_types.csv",header=TRUE)
plasma_genetype <- read.csv("C:/Users/Tao/Desktop/Candidates/STAD/相关数据/saliva GC biomarker/plasma_gene_types.csv",header=TRUE)

ggplot(plasma_genetype,aes(x="",fill=Types))+
  geom_bar(stat="Count",width=0.5,position='stack')+
  coord_polar("y", start=0)+
  xlab("")+
  ylab("")+
  geom_text(stat="count",aes(label = scales::percent(..count../13729)), size=4, position=position_stack(vjust = 0.5))+
  theme_bw()+
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype="solid",,size=0.7,colour = "white"))
