#chrM boxplot
{
test <- read.csv("/Users/yuhuan/Desktop/genome.txt",sep = "\t", header = TRUE, row.names = 1)
test <- as.data.frame(cpm(test))
plot <- as.data.frame(t(test[which(rownames(test)=="chrM"),]))
plot$group <- as.character(lapply(strsplit(rownames(plot),".",fixed=TRUE),function(x) x[1]))

ggplot(plot,aes(x=group,y=chrM))+geom_boxplot()
  
my_comparisons <- list(c("NC","CRC"),c("NC","STAD"))
ggplot(plot,aes(x=group,y=chrM))+
  geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
  #geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
  scale_fill_brewer(palette="Blues") +
  theme_bw()+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(size=1, colour = "black"),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
    axis.text.y = element_text(face="bold",  color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  stat_compare_means(comparisons = my_comparisons,
                       method = "wilcox.test",
                       method.args = list(alternative = "two.sided",paired = TRUE),
                       label = "p.signif")
}

#Expression candidates
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
Differential_result <- read.csv("./Expression/output/Expression_STADvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

Differential_result$ensembl_gene_id <- as.character(lapply(strsplit(rownames(Differential_result),"\\."),function(x) x[1]))

gene_names <- getBM(attributes=c("ensembl_gene_id", "external_gene_name","chromosome_name"),
                        filters = "ensembl_gene_id",
                        values=Differential_result$ensembl_gene_id, mart= mart,useCache = FALSE)
Differential_result$Rownames <- rownames(Differential_result)
Differential_result_new <- left_join(Differential_result,gene_names,by = c("ensembl_gene_id"="ensembl_gene_id"))

TPM_Expression <- read.csv("./Expression/matrix/20211113_STAD_vs_NC_TPM.txt",header = T, row.names = 1,sep = "\t")

gini_Expression_STAD <- gini(t(TPM_Expression[,grep("STAD",colnames(TPM_Expression))]))
gini_Expression_STAD <- as.data.frame(gini_Expression_STAD)
colnames(gini_Expression_STAD) <- c("gini_index_STAD")
gini_Expression_NC <- gini(t(TPM_Expression[,grep("NC",colnames(TPM_Expression))]))
gini_Expression_NC <- as.data.frame(gini_Expression_NC)
colnames(gini_Expression_NC) <- c("gini_index_NC")
gini_Expression <- cbind(gini_Expression_NC,gini_Expression_STAD)

mean_Expression_STAD <- rowMeans(TPM_Expression[,grep("STAD",colnames(TPM_Expression))])
mean_Expression_STAD <- as.data.frame(mean_Expression_STAD)
colnames(mean_Expression_STAD) <- c("mean_TPM_STAD")
mean_Expression_NC <- rowMeans(TPM_Expression[,grep("NC",colnames(TPM_Expression))])
mean_Expression_NC <- as.data.frame(mean_Expression_NC)
colnames(mean_Expression_NC) <- c("mean_TPM_NC")
mean_Expression <- cbind(mean_Expression_NC,mean_Expression_STAD)

gini_TPM <- cbind(gini_Expression,mean_Expression,TPM_Expression)

#rownames(Differential_result_new) <- Differential_result_new$Rownames
#Differential_result_new$rownames <- rownames(Differential_result_new)
gini_TPM$rownames <- rownames(gini_TPM)

Differential_result_new2 <- left_join(Differential_result_new,gini_TPM,by = c("Rownames"="rownames"))
Differential_result_new2 <- Differential_result_new2[,-which(colnames(Differential_result_new2)=="rownames")]

rownames(Differential_result_new2) <- rownames(Differential_result_new)
write.table(Differential_result_new2,"Expression_STADvsNC_edger_exact_Name_Chr.txt",sep = "\t", quote = FALSE)
}

#Expression candidates
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
  Differential_result <- read.csv("./Medip/output/STADvsNC-passQCstringent-edger_exact.txt",header = T, row.names = 1,sep = "\t")
  rownames(Differential_result) <- gsub("promoter_","",rownames(Differential_result))
  library(biomaRt)
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  
  Differential_result$ensembl_gene_id <- as.character(lapply(strsplit(rownames(Differential_result),"\\."),function(x) x[1]))
  
  gene_names <- getBM(attributes=c("ensembl_gene_id", "external_gene_name","chromosome_name"),
                      filters = "ensembl_gene_id",
                      values=Differential_result$ensembl_gene_id, mart= mart,useCache = FALSE)
  Differential_result$Rownames <- rownames(Differential_result)
  Differential_result_new <- left_join(Differential_result,gene_names,by = c("ensembl_gene_id"="ensembl_gene_id"))
  
  TPM_Expression <- read.csv("./Medip/matrix/CPMtmm_matrix_promoter.txt",header = T, row.names = 1,sep = "\t")
  
  gini_Expression_STAD <- gini(t(TPM_Expression[,grep("STAD",colnames(TPM_Expression))]))
  gini_Expression_STAD <- as.data.frame(gini_Expression_STAD)
  colnames(gini_Expression_STAD) <- c("gini_index_STAD")
  gini_Expression_NC <- gini(t(TPM_Expression[,grep("NC",colnames(TPM_Expression))]))
  gini_Expression_NC <- as.data.frame(gini_Expression_NC)
  colnames(gini_Expression_NC) <- c("gini_index_NC")
  gini_Expression <- cbind(gini_Expression_NC,gini_Expression_STAD)
  
  mean_Expression_STAD <- rowMeans(TPM_Expression[,grep("STAD",colnames(TPM_Expression))])
  mean_Expression_STAD <- as.data.frame(mean_Expression_STAD)
  colnames(mean_Expression_STAD) <- c("mean_CPM_STAD")
  mean_Expression_NC <- rowMeans(TPM_Expression[,grep("NC",colnames(TPM_Expression))])
  mean_Expression_NC <- as.data.frame(mean_Expression_NC)
  colnames(mean_Expression_NC) <- c("mean_CPM_NC")
  mean_Expression <- cbind(mean_Expression_NC,mean_Expression_STAD)
  
  gini_TPM <- cbind(gini_Expression,mean_Expression,TPM_Expression)
  
  #rownames(Differential_result_new) <- Differential_result_new$Rownames
  #Differential_result_new$rownames <- rownames(Differential_result_new)
  gini_TPM$rownames <- gsub("promoter_","",rownames(gini_TPM))
  
  Differential_result_new2 <- left_join(Differential_result_new,gini_TPM,by = c("Rownames"="rownames"))
  Differential_result_new2 <- Differential_result_new2[,-which(colnames(Differential_result_new2)=="Rownames")]
  
  Differential_result_new2$Rownames <- Differential_result_new$Rownames
  write.table(Differential_result_new2,"MeDIP_STADvsNC_edger_exact_Name_Chr.txt",sep = "\t", quote = FALSE)
}

#miRNA candidates
{
Differential_result_miRNA <- read.csv("./miRNA/output/miRNA_STADvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")

Differential_result_miRNA$ensembl_transcript_id <- as.character(lapply(strsplit(rownames(Differential_result_miRNA),"\\."),function(x) x[1]))

gene_names_miRNA <- getBM(attributes=c("ensembl_transcript_id", "external_gene_name","chromosome_name"),
                    filters = "ensembl_transcript_id",
                    values=Differential_result_miRNA$ensembl_transcript_id, mart= mart,useCache = FALSE)

Differential_result_miRNA_new <- left_join(Differential_result_miRNA,gene_names_miRNA,by = c("ensembl_transcript_id"="ensembl_transcript_id"))

rownames(Differential_result_miRNA_new) <- rownames(Differential_result_miRNA)

Count_miRNA <- read.csv("./miRNA/matrix/miRNA_for_STADvsNC.txt",header = T, row.names = 1,sep = "\t")

CPM_miRNA <- cpm(Count_miRNA)

gini_miRNA_STAD <- gini(t(CPM_miRNA[,grep("STAD",colnames(CPM_miRNA))]))
gini_miRNA_STAD <- as.data.frame(gini_miRNA_STAD)
colnames(gini_miRNA_STAD) <- c("gini_index_STAD")
gini_miRNA_NC <- gini(t(CPM_miRNA[,grep("NC",colnames(CPM_miRNA))]))
gini_miRNA_NC <- as.data.frame(gini_miRNA_NC)
colnames(gini_miRNA_NC) <- c("gini_index_NC")
gini_miRNA <- cbind(gini_miRNA_NC,gini_miRNA_STAD)

mean_miRNA_STAD <- rowMeans(CPM_miRNA[,grep("STAD",colnames(CPM_miRNA))])
mean_miRNA_STAD <- as.data.frame(mean_miRNA_STAD)
colnames(mean_miRNA_STAD) <- c("mean_CPM_STAD")
mean_miRNA_NC <- rowMeans(CPM_miRNA[,grep("NC",colnames(CPM_miRNA))])
mean_miRNA_NC <- as.data.frame(mean_miRNA_NC)
colnames(mean_miRNA_NC) <- c("mean_CPM_NC")
mean_miRNA <- cbind(mean_miRNA_NC,mean_miRNA_STAD)

gini_CPM_miRNA <- cbind(gini_miRNA,mean_miRNA,CPM_miRNA)

rownames(Differential_result_miRNA_new) <- rownames(Differential_result_miRNA)
Differential_result_miRNA_new$rownames <- rownames(Differential_result_miRNA_new)
gini_CPM_miRNA$rownames <- rownames(gini_CPM_miRNA)

Differential_result_miRNA_new2 <- left_join(Differential_result_miRNA_new,gini_CPM_miRNA,by = c("rownames"="rownames"))
Differential_result_miRNA_new2 <- Differential_result_miRNA_new2[,-which(colnames(Differential_result_miRNA_new2)=="rownames")]

rownames(Differential_result_miRNA_new2) <- rownames(Differential_result_miRNA_new)

write.table(Differential_result_miRNA_new2,"miRNA_STADvsNC_edger_exact_Name_Chr.txt",sep = "\t", quote = FALSE)
}

#lncRNA-GC1 and immune gene boxplot
{
test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20210722_multiomics/Expression_20211113/Expression_ML.txt",sep = "\t", header = TRUE, row.names = 1)
#plot <- as.data.frame(t(test[which(rownames(test)=="lncRNA-GC1|2145"),]))
#test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20210722_multiomics/Altpromoter/Altpromoter_ML.txt",sep = "\t", header = TRUE, row.names = 1)
plot <- as.data.frame(t(test[grep("ENSG00000198793",rownames(test)),]))


#wangtantan HPK1, MAP4K1, ENSG00000104814
#mTOR: ENSG00000198793

#Upregulated inhibitory receptor
#ENSG00000120217: CD274(PD-L1)
#ENSG00000103855: CD276(B7-H3)
#ENSG00000134258: VTCN1(B7-H4)
#ENSG00000114455: HHLA2(B7-H5)
#ENSG00000188389: PDCD1(PD1)
#ENSG00000163599: CTLA4
#ENSG00000089692: LAG3
#ENSG00000135077: TIM3
#ENSG00000079385: CEACAM1(TIM3 ligand)
#ENSG00000181847: TIGIT

#Downregulated immune activation genes
#ENSG00000153563: CD8A
#ENSG00000172116: CD8B
#ENSG00000010610: CD4
#ENSG00000115085: ZAP70
#ENSG00000113263: ITK
#ENSG00000111537: IFNG
#ENSG00000074966: RLK
#ENSG00000170345: FOS
plot$group <- as.character(lapply(strsplit(rownames(plot),".",fixed=TRUE),function(x) x[1]))

plot$group <- gsub("NC","HD",plot$group)

my_comparisons <- list(c("HD","CRC"),c("HD","STAD"))
plot$group <- factor(plot$group,levels=c("HD","CRC","STAD"))

plot$group <- gsub("CRC","GIC",plot$group)
plot$group <- gsub("STAD","GIC",plot$group)
my_comparisons <- list(c("HD","GIC"))
plot$group <- factor(plot$group,levels=c("HD","GIC"))
forplot <- plot[,c(1,ncol(plot))]
colnames(forplot) <- c("value","group")
ggplot(forplot[-which(rownames(forplot)=="CRC.PKU.29.pico"),],aes(x=group,y=value,fill = group))+
#ggplot(forplot,aes(x=group,y=value,fill = group))+
  geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
  geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
  scale_fill_manual(values=c("#87CEEB","#FCB514"))+
  #scale_fill_brewer(palette="Blues") +
  #ylim(0,25)+
  theme_bw()+
  xlab("")+
  ylab("mTOR")+
  #ylab(colnames(plot)[1])+
  theme(#legend.position="right",
    legend.position="none",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(size=1, colour = "black"),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
    axis.text.y = element_text(face="bold",  color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "greater",paired = TRUE),
                     label = "p.signif",
                     size = 10,
                     vjust = 0.5)
}

#ML plot (for shaozhen qualify)
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20211114_yunfan/")
plot <- read.csv("combination-performance.txt",header = TRUE,sep = "\t")

plot <- plot[(plot$metric=="AUROC")&(plot$cancer=="STAD-NC"),]

#plot ROC-curve
{
  library(pROC)
  #report ROC
  files <- grep("-",dir("cross-validation-combination/STAD-NC/7/4/"),value = TRUE)
  files <- grep("expression",files,value = TRUE)
  files <- grep("MeDIP",files,value = TRUE)
  files <- grep("nucleasome",files,value = TRUE)
  files <- grep("CNV",files,value = TRUE)
  i=1
  while(i<=length(files)){
    predicted <- read.csv(paste0("cross-validation-combination/STAD-NC/7/4/",files[i]),sep = "\t")
    roc.curve <- roc(predicted$y_true,predicted$y_pred,levels=c(0,1),direction=c("<"))
    ci.auc(roc.curve,conf.level = 0.95)
    record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
    print(paste0(files[i],":",roc.curve$auc))
    i=i+1
  }

ML_samples <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/ML_samples.csv",header = TRUE)
predicted_all <- read.csv("cross-validation-combination/STAD-NC/5/4/MeDIP-expression-nucleasome-CNV.txt",sep = "\t")

predicted_meta <- left_join(predicted_all,ML_samples,by = c("sample_id"="sample_id"))

#alteration
{
  predicted_expression <- read.csv("cross-validation-combination/STAD-NC/5/1/expression.txt",sep = "\t")
  predicted_medip <- read.csv("cross-validation-combination/STAD-NC/5/1/MeDIP.txt",sep = "\t")
  predicted_nucleasome <- read.csv("cross-validation-combination/STAD-NC/5/1/nucleasome.txt",sep = "\t")
  predicted_cnv <- read.csv("cross-validation-combination/STAD-NC/5/1/CNV.txt",sep = "\t")
  
  RNA_color <- c("Alternative promoter"="#E9C2A6",
                 "RNA Expression"="#A5435C",
                 "RNA Splicing"="#C5E3BF",
                 "Alternative polyadenyltion"="#003F87",
                 "Chimeric RNA"="#FF3D0D",
                 "RNA editing"="#324F17",
                 "Allele specific expression"="#800080",
                 "RNA SNP"="#333333",
                 "DNA copy number"="#FFE234",
                 "DNA nucleosome occupancy"="#44A5F9",
                 "DNA methylation"="#C7C2C2",
                 "Combined"="white")
  
  par(pty="s")
  roc(predicted_all$y_true,predicted_all$y_pred,plot=TRUE,legacy.axes=TRUE,
      xlab="False Positive Percentage",ylab="True Positive Percentage",
      col="purple", lwd=4,print.auc=TRUE,print.auc.y= 0.30,print.auc.x= 0.45,cex.axis = 1.5, cex.lab = 2)
  plot.roc(predicted_expression$y_true,predicted_expression$y_pred,
           col="#A5435C", lwd=4, print.auc=TRUE, print.auc.y=0.25, print.auc.x= 0.45, add=TRUE)
  plot.roc(predicted_medip$y_true,predicted_medip$y_pred,
           col="#C7C2C2", lwd=4, print.auc=TRUE, print.auc.y=0.20, print.auc.x= 0.45, add=TRUE)
  
  plot.roc(predicted_nucleasome$y_true,predicted_nucleasome$y_pred,
           col="#44A5F9", lwd=4, print.auc=TRUE, print.auc.y=0.15, print.auc.x= 0.45, add=TRUE)
  plot.roc(predicted_cnv$y_true,predicted_cnv$y_pred,
           col="#FFE234", lwd=4, print.auc=TRUE, print.auc.y=0.05, print.auc.x= 0.45, add=TRUE)
  legend("bottomright",legend=c("4-Combined: 0.910","Nucleasome occupancy: 0.854","RNA expression:0.776","DNA methylation:0.668","DNA copy number: 0.546"),
         col=c("purple","#44A5F9","#A5435C","#C7C2C2","#FFE234"),lwd=4,
         cex = 1.2)
}

#stage
{
  predicted_1 <- predicted_meta[(predicted_meta$Stage_numeric==1)|(predicted_meta$Stage_numeric==0),c(1,2,3)]
  predicted_2 <- predicted_meta[(predicted_meta$Stage_numeric==2)|(predicted_meta$Stage_numeric==0),c(1,2,3)]
  predicted_3 <- predicted_meta[(predicted_meta$Stage_numeric==3)|(predicted_meta$Stage_numeric==0),c(1,2,3)]
  predicted_4 <- predicted_meta[(predicted_meta$Stage_numeric==4)|(predicted_meta$Stage_numeric==0),c(1,2,3)]
  
  par(pty="s")
  roc(predicted_1$y_true,predicted_1$y_pred,plot=TRUE,legacy.axes=TRUE,
      xlab="False Positive Percentage",ylab="True Positive Percentage",
      col="red", lwd=4,print.auc=TRUE,print.auc.y= 0.20,print.auc.x= 0.35,cex.axis = 1.5, cex.lab = 2)
  plot.roc(predicted_2$y_true,predicted_2$y_pred,
           col="#EE9A00", lwd=4, print.auc=TRUE, print.auc.y=0.15, print.auc.x= 0.35, add=TRUE)
  plot.roc(predicted_3$y_true,predicted_3$y_pred,
           col="grey", lwd=4, print.auc=TRUE, print.auc.y=0.10, print.auc.x= 0.35, add=TRUE)
  plot.roc(predicted_4$y_true,predicted_4$y_pred,
           col="black", lwd=4, print.auc=TRUE, print.auc.y=0.05, print.auc.x= 0.35, add=TRUE)
  legend("bottomright",legend=c("Stage I: 0.952","Stage II: 0.875","Stage III: 0.922","Stage IV: 0.838"),
         col=c("red","#EE9A00","grey","black"),lwd=4,cex=1.6)
}

#gender
{
  predicted_male <- predicted_meta[predicted_meta$Gender=="M",c(1,2,3)]
  predicted_female <- predicted_meta[(predicted_meta$Gender=="F"),c(1,2,3)]
  
  par(pty="s")
  roc(predicted_all$y_true,predicted_all$y_pred,plot=TRUE,legacy.axes=TRUE,
      xlab="False Positive Percentage",ylab="True Positive Percentage",
      col="purple", lwd=4,print.auc=TRUE,print.auc.y= 0.20,print.auc.x= 0.35,cex.axis = 1.5, cex.lab = 2)
  plot.roc(predicted_male$y_true,predicted_male$y_pred,
           col="dark blue", lwd=4, print.auc=TRUE, print.auc.y=0.15, print.auc.x= 0.35, add=TRUE)
  plot.roc(predicted_female$y_true,predicted_female$y_pred,
           col="pink", lwd=4, print.auc=TRUE, print.auc.y=0.10, print.auc.x= 0.35, add=TRUE)
  legend("bottomright",legend=c("All gender: 0.914","Male: 0.994","Female: 0.754"),
         col=c("purple","dark blue","pink"),lwd=4,cex=1.6)
}

#Age
{
  predicted_older <- predicted_meta[predicted_meta$Age>=60,c(1,2,3)]
  predicted_younger <- predicted_meta[(predicted_meta$Age<60),c(1,2,3)]
  
  par(pty="s")
  roc(predicted_all$y_true,predicted_all$y_pred,plot=TRUE,legacy.axes=TRUE,
      xlab="False Positive Percentage",ylab="True Positive Percentage",
      col="purple", lwd=4,print.auc=TRUE,print.auc.y= 0.20,print.auc.x= 0.35,cex.axis = 1.5, cex.lab = 2)
  plot.roc(predicted_older$y_true,predicted_older$y_pred,
           col="brown", lwd=4, print.auc=TRUE, print.auc.y=0.15, print.auc.x= 0.35, add=TRUE)
  plot.roc(predicted_younger$y_true,predicted_younger$y_pred,
           col="green", lwd=4, print.auc=TRUE, print.auc.y=0.10, print.auc.x= 0.35, add=TRUE)
  legend("bottomright",legend=c("All: 0.914","Older(>=60): 0.861","Younger(<60): 0.940"),
         col=c("purple","brown","green"),lwd=4,cex=1.6)
}
}

#plot features of a specific combination(50 features)
{
  features_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20211114_yunfan/features"
  group <- "cancer-NC"
  alterations <- c("MeDIP","alt-promoter","editing","nucleasome","SNP","chimera")
  number <- length(alterations)
  seed <- "5"
  
  #summarize features for a specific combination
  i=1
  features <- {}
  while(i<=length(alterations)){
    features_LOOCV <- dir(paste0(features_dir,"/",group,"/",alterations[i],"/",seed))
    j=1
    single_features <- {}
    while(j<=length(features_LOOCV)){
      single_alteration_features <- read.csv(paste0(features_dir,"/",group,"/",alterations[i],"/",seed,"/",features_LOOCV[j]),sep = "\t",header = FALSE)
      single_alteration_features <- head(single_alteration_features,floor(50/number))
      single_alteration_features$type <- alterations[i]
      single_features <- rbind(single_features,single_alteration_features)
      j=j+1
    }
    features <- rbind(features,single_features)
    i=i+1
  }
  
  
  #
  features$count <- 1
  colnames(features) <- c("alteration","importance","type","count")
  features$alteration <- paste(features$type,features$alteration,sep = "|")
  feature_importance <- aggregate(importance ~ alteration, features, mean)
  feature_reccurence <- aggregate(count ~ alteration, features, sum)
  
  feature_plot <- left_join(feature_reccurence,feature_importance, by = c("alteration"="alteration"))
  
  final_plot <- feature_plot[feature_plot$count>=(length(features_LOOCV))*0.5,]
  
  final_plot <- head(final_plot[order(final_plot$importance,decreasing = TRUE),],20)
  #write.csv(final_plot,"demo_feature.csv")
  final_plot <- read.csv("demo_feature.csv")
  final_plot$name <- factor(final_plot$name,levels = final_plot$name)
  final_plot <- head(final_plot[order(final_plot$importance,decreasing = TRUE),],15)
  ggplot(final_plot,aes(x=importance,y=name,fill=-log10(pvalue)))+
    geom_bar(stat='identity')+
    scale_colour_gradient2(low="#162252",high="#3B5E2B",mid="#008B00",midpoint = 0)+
    #geom_vline(xintercept = c(6.5,10.5,16.5),linetype = "dashed",color="black")+
    #geom_vline(xintercept = 25.5,linetype = "dashed",color="grey")+
    labs(color="Importance",
         #size=expression("Reccurence"),
         x="")+
    theme_bw()+
    theme(strip.text = element_text(size = rel(1.3),face="bold",colour = "black"),
          axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
          axis.text.x = element_text(size=rel(1.3),face="bold",colour = "black",angle = 45,vjust=1,hjust=1),
          axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
          axis.title.y = element_blank())
  
}

#sparseness
{
  path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20211105_yunfan/ML_matrix"
  
  alterations <- gsub("_ML.txt","",dir(path),fixed=TRUE)
  alteration_sparseness <- data.frame(matrix(ncol=3))
  colnames(alteration_sparseness) <- c("alteration","sparseness","alteration_number")
  i=1
  while(i<=length(alterations)){
    matrix <- read.csv(paste0(path,"/",alterations[i],"_ML.txt"),sep = "\t",header = TRUE)
    matrix[is.na(matrix)] <- 0
    total <- (nrow(matrix)*ncol(matrix))
    zero <- sum(matrix == 0)
    sparseness <- zero/total
    alteration_sparseness[i,"alteration"] <- alterations[i]
    alteration_sparseness[i,"sparseness"] <- sparseness
    alteration_sparseness[i,"alteration_number"] <- nrow(matrix)
    i=i+1
  }
  alteration_sparseness$`1-sparseness` <- 1-alteration_sparseness$sparseness
  write.csv(alteration_sparseness,"alteration_sparseness.csv")
}

SF1 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/Figures/V3/SF1.csv",header =TRUE)
SF1$mixed_sample <- factor(SF1$mixed_sample,levels = rev(as.character(SF1$mixed_sample)))
ggplot(SF1,aes(x=number,y=mixed_sample,fill=type))+
  geom_bar(stat='identity', color = "black")+
  #scale_fill_gradient2(low="#B0E2FF",high="#162252",mid="#26466D",midpoint = 5,na.value = "white")+
  #scale_x_continuous(expand = c(0,0.0001))+
  #geom_vline(xintercept = c(6.5,10.5,16.5),linetype = "dashed",color="black")+
  #geom_vline(xintercept = 25.5,linetype = "dashed",color="grey")+
  labs(fill="Sample type",
       x="")+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = rel(1.3),face="bold",colour = "black"),
        axis.text.y = element_text(size = rel(1),colour = "black"),
        axis.text.x = element_text(size=rel(1),colour = "black"),
        axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
        axis.title.y = element_blank())
}

