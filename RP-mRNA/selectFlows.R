#### environment preparation
{
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(tidyverse)
  library(reshape)
  library(extrafont)
  library(edgeR)
  library(DESeq2)
  fonts()
  library(progress)
  library(pheatmap)
  library(clusterProfiler) ## enrichment
  library(biomaRt) ## ID convert
  library(org.Hs.eg.db) ## annotation for human
  library(enrichplot) ##gseaplot
  library(GO.db)
}

#pathway preperation
{
  #kegg
  hsa_kegg <- clusterProfiler::download_KEGG("hsa")
  names(hsa_kegg)
  PATH2NAME <- hsa_kegg$KEGGPATHID2NAME
  PATH2ID <- hsa_kegg$KEGGPATHID2EXTID
  PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by="from")
  colnames(PATH_ID_NAME) <- c("KEGGID", "ENTREZID", "DESCRPTION")
  #ID conversion
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  entrezgene <- PATH_ID_NAME$ENTREZID
  ensembl_gene_id<- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                          filters = "entrezgene_id",
                          values=entrezgene, mart= mart,useCache = FALSE)
  PATH_ID_NAME <- merge(PATH_ID_NAME, ensembl_gene_id, by.x= "ENTREZID",by.y= "entrezgene_id")
  ensg <- PATH_ID_NAME$ensembl_gene_id
  gene_symbol <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),
                       filters = "ensembl_gene_id",
                       values=ensg, mart= mart,useCache = FALSE)
  PATH_ID_NAME <- merge(PATH_ID_NAME, gene_symbol, by.x= "ensembl_gene_id",by.y= "ensembl_gene_id")
  PATH_ID_NAME <- unique(PATH_ID_NAME)
  
  #read from local
  PATH_ID_NAME <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/PATH_ID_NAME.csv",header = TRUE,row.names=1)
  PATH_ID_NAME <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/PATH_ID_NAME_modified.csv",header = TRUE,row.names=1)
  
  Ribosome <- subset(PATH_ID_NAME,DESCRPTION=="Ribosome")
  Ribosome_ensg <- unique(Ribosome$ensembl_gene_id)
  Ribosome_symbol <- unique(Ribosome$hgnc_symbol)
  
  #GO
  GO <- getBM(attributes=c("ensembl_gene_id","go_id"),
                       filters = "ensembl_gene_id",
                       values=ensg, mart= mart,useCache = FALSE)
  columns(GO.db)
  keytypes(GO.db)
  keys <- head(keys(GO.db))
  keys
  GO_ID_NAME <- select(GO.db, keys = keys(GO.db),columns = c("GOID","TERM","ONTOLOGY"))
  GO_ID_NAME <- merge(GO_ID_NAME, GO, by.x= "GOID",by.y= "go_id")
}


#### differential analysis
### DESeq2
{
  #DESeq2的步骤就3步：1.构建dds矩阵 2.标准化 3.差异分析
  
  #~在R里面用于构建公式，~左边为因变量，右边为自变量
  mat <- read.table("./matrix/counts-G_filtered.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据（txt类型）
  des <- read.csv("./DE/des_CRCvsSTAD.csv", header = TRUE, check.names=FALSE, sep=',') #读入数据（csv类型）
  samples <- des$samples   #取出来每一列
  group <- des$group
  batch <- des$batch
  des$group <- as.factor(des$group)
  des$batch <- as.factor(des$batch)
  dds <- DESeqDataSetFromMatrix(countData = mat,    #构建dds矩阵：countData需要样本信息表，colData需要样品信息矩阵
                                colData = as.data.frame(des),  #样品信息矩阵(格式要求是data.frame)，样品名称不变样品处理情况需要是factor
                                design = ~group+batch)  #差异表达矩阵：差异比较矩阵就是告诉差异分析函数是要从要分析哪些变量间的差异，简单说就是说明哪些是对照哪些是处理                
  
  #加个batch是啥意思：和group一样都是变量
  dds <- DESeq(dds,test="Wald")               #对dds进行标准化，wald检测
  res <- results(dds, contrast=c('group', 'positive', 'negative'))
  summary(res) #查看一下分析的摘要
  
  write.table(as.data.frame(res), "deseq2_CRCvsSTAD.csv", sep=',', quote=FALSE, row.names=TRUE)   #输出结果
}    

### EdgeR_glmlrt
{
  mat <- read.table("./matrix/counts-G_filtered.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
  des <- read.csv("./DE/des_CRCvsSTAD.csv", header = TRUE, check.names=FALSE, sep=',')  #读入数据
  samples <- des$samples
  group <- des$group
  batch <- des$batch
  y <- DGEList(counts=mat, samples=samples, group=group)     #规定输入格式，必须包含这几个参数
  y <- calcNormFactors(y, method="TMM")                      #归一化处理
  design <- model.matrix(~group+batch)   #group是区别分组的factor的不同样本，lane是不同样本的lane道，意思是group和batch 是自变量进行分析  
  y <- estimateDisp(y, design)          #广义线性模型计算离散度（common&trended&tagwise）
  fit <- glmFit(y, design)              #差异表达分析函数                     
  test <- glmLRT(fit, coef=2)           #差异表达分析函数 
  res <- topTags(test, n=nrow(mat), sort.by='none')     #输出计算差异基因的结果
  res <- cbind(res$table, baseMean=2^(res$table$logCPM)) #输出差异基因结果加上baseMean一栏
  mapped_names <- colnames(res)                   #后面都是整理+输出
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  write.table(res, "edger_glmlrt_CRCvsSTAD.csv", sep=',', quote=FALSE, row.names=TRUE) #输出文件，规定名字
}


#### correclation between deseq2 and edger_glmlrt
{
  deseq <- read.csv("deseq2_CRCvsSTAD.txt", sep = "\t", header =TRUE, row.names = 1)  #输入文件1
  egder <- read.csv("edger_glmlrt_CRCvsSTAD.txt", sep = "\t", header =TRUE, row.names = 1) #输入文件2
  
  head(egder)
  
  ggplot(deseq,aes(deseq$log2FoldChange)) +
    geom_histogram()              #看一下数据情况，是上调的多下调的多？以及其他信息，主要起预览作用
  
  FC={}
  FC$egder <- egder$log2FoldChange    
  FC$deseq <- deseq$log2FoldChange
  FC <- as.data.frame(FC)
  rownames(FC) <- row.names(deseq)                   #以上为构建FC表格进行分析两者相关性
  
  r <- cor(FC$egder,FC$deseq,method="pearson")       #相关性分析（皮尔森）也就是后面标注的R相关系数
  ggplot(data=FC, aes(x=egder, y=deseq)) +           #后面就是画散点图分析来相关性
    geom_point(color="#d7191c") +
    #geom_smooth(method="lm",color="#1a9641") +
    geom_text(aes(x=2, y=-3,label=paste("R","=",signif(r,3),seq="")),color="#fdae61",size=12)+ 
    #这一步有标识R的一步：用geom_text添加，设置包括：位置、和拼接内容，颜色，字体大小
    theme_bw()+
    xlab("Edger")+
    ylab("DESeq2")+
    theme(
      axis.title = element_text(face="bold", color="black", size=24),
      axis.text = element_text(face="bold",  color="black", size=24)
    )
}


#### figure Differential expression between NC_PKU and othersNC
{
  volcano <- read.csv("deseq2_CRCvsNC_PKU_decontamed.txt", sep = '\t', header = T, row.names = 1, fileEncoding = 'utf-8', fill = T)     #导入文件
  volcano$threshold <- as.factor(ifelse(volcano$padj < 0.05 & abs(volcano$log2FoldChange) >=1,ifelse(volcano$log2FoldChange > 1 ,'Up','Down'),'Not'))   #设置标记标签(factor类型）
  write.csv(volcano,"CRCvsSTAD.csv")    #输出文件
  head(volcano)
  #volcano <- volcano[complete.cases(volcano), ]
  #head(volcano)
  #rownames(volcano) <- volcano[,1]
  volcano_plot <- ggplot(data=volcano, aes(x=log2FoldChange, y =-log10(pvalue), colour=threshold)) + #颜色根据threshold分（是factor）
    scale_color_manual(values=c("blue", "grey","red"))+ #设置颜色
    geom_point() + #画散点
    #xlim(c(-6,8.5)) +
    #ylim(c(0,5)) +
    theme_bw(base_size = 12, base_family = "Arial") + #装饰
    geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6,aes(0,4))+ #阀线
    geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6,aes(-3,6))+ #阀线
    theme(#legend.position="right", #装饰（下面都是）
      legend.position="right",
      panel.grid=element_blank(),
      #legend.title = element_blank(),
      #legend.text= element_text(face="bold", color="black",family = "Times", size=8),
      plot.title = element_text(hjust = 0.5,size=30,face="bold"),
      axis.text.x = element_text(face="bold",family = "Arial",color="black", size=30),
      axis.text.y = element_text(face="bold",family = "Arial",color="black", size=30),
      axis.title.x = element_text(face="bold",family = "Arial",color="black", size=36),
      axis.title.y = element_text(face="bold",family = "Arial",color="black", size=36))+
    labs(x="log2(fold change)",y="-log10 (p.adj)",title="Different Microbe Abundance", face="bold")
  volcano_plot   
  
}


#### boxplot to illustrate batch effect (raw abundance normalized to relative abundance)

## read in files
{
  #raw to illustrate batch
  counts <- read.csv("./matrix/counts-G_filtered_norm.txt",sep="\t",header = T, row.names = 1)
  
  #check distribution
  head(grep("Cyprinivirus",rownames(counts),value = T),10)  #value如果是F返回的是索引，T返回的是本身
  distribution <- as.data.frame(t(counts))                  #倒转counts          
  ggplot(distribution,aes(distribution$`1279|Staphylococcus`)) +    #画图看看这个菌在各个微生物中表现如何
    geom_histogram()
}

##manually vlookup genome size

## input candidates
{
  genome <- read.csv("CRCvsNC_PKU_genome.csv",header = T)
  #candidate <- rownames(counts)
  candidate <- genome$candidate       #输入candidates的基因名字
  genome_size <- genome$size          #输入candidates的基因组大小（准备进行标准化）
  logFC <- genome$logFC
}

####  microbe selection workflow
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/06.selectbywilcox_microbial/")
  counts <- read.csv("counts-G_norm_decontamed.txt",sep="\t",header = T, row.names = 1)
  positive <- c("CRC")
  negative <- c("NC_PKU","STAD","HCC","LUAD","ESCA")
  candidate <- rownames(counts)
  n=1
  #matrix_log2counts = log2(counts+1)
  while(n<=length(negative)){
    #logFC <- apply(matrix_log2counts[,grep(positive[1],colnames(matrix_log2counts),fixed=TRUE)], 1, mean)-apply(matrix_log2counts[,grep(negative[n],colnames(matrix_log2counts),fixed=TRUE)], 1, mean)
    #logFC <- as.data.frame(logFC)
    j=1
    #while(j<=length(candidate)){
    #  print(grep(candidate[j],rownames(counts),fixed=TRUE))
    #  j=j+1
    #}
    test={}
    while(j<=length(candidate)){
      target <- candidate[j]
      type <- c(positive[1],negative[n])
      i=1
      gene={}
      while(i<=length(type)){
        temp <- counts[grep(target,rownames(counts),fixed=TRUE)[1],grep(type[i],colnames(counts),fixed=TRUE)]
        temp.t <- t(temp)
        colnames(temp.t) <- c("counts")
        temp.df <- as.data.frame(temp.t)
      #remove outlier
        #k=1
        #OutVals <- boxplot(temp.t,plot=FALSE)$out
        #OutVals <- unique(OutVals)
        #OutVals <- sort(OutVals,decreasing = T)
        #while(k<=length(OutVals)){
        #  # return integer(0), temp.df <- temp.df
        #  if(length(grep(OutVals[k],temp.df$counts))==0){
        #   temp.df <- temp.df
        #  } else if(OutVals[k]==0) {
        #    temp.df <- temp.df
        #  }
        #  else {
        #    temp.df <- as.data.frame(temp.df[-grep(OutVals[k],temp.df$counts,fixed=TRUE),])
        #    colnames(temp.df) <- "counts"
        #    #print(length(temp.df$counts))
        #  }
        #  k=k+1
        #}
        temp.df$CancerType<-rep(type[i],times=length(temp.df$counts))
        temp.df$microbe<-rep(target,times=length(temp.df$counts))
        gene <- rbind(gene,temp.df)
        
        i=i+1
      }
      trend <- mean(subset(gene,CancerType==positive[1])$counts)-mean(subset(gene,CancerType==negative[n])$counts)
      
      #wilcox.test       
      if(trend>0){
        #p <- p+stat_compare_means(comparisons = my_comparisons,
        #                          method = "wilcox.test",
        #                          method.args = list(alternative = "greater"),
        #                          label = "p.signif")+labs(x="",y="relative abundance",title=paste0(target,"_wilcox.test.greater"), face="bold")
        #ggsave(p,filename = paste0(target,"_greater.pdf"),path = "./")
        t <- wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="greater" )
      } else if(trend<0) {
        #p <- p+stat_compare_means(comparisons = my_comparisons,
        #                          method = "wilcox.test",
        #                          method.args = list(alternative = "less"),
        #                          label = "p.signif")+labs(x="",y="relative abundance",title=paste0(target,"_wilcox.test.less"), face="bold")
        #ggsave(p,filename = paste0(target,"_less.pdf"),path="./")
        t <-wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="less")
      } else {
        t <-wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="two.sided")
      }
      t$microbe <- rep(target,times=length(t$p.value)) 
      t$positive_vs_negative <- rep(paste0(positive[1],"_vs_",negative[n]),times=length(t$p.value))
      test <- rbind(test,t)
      j=j+1
    }
    t.df={}
    test <- as.data.frame(test)
    t.df$microbe <- unlist(test$microbe)
    t.df$positive_vs_negative <- unlist(test$positive_vs_negative)
    t.df$p.value <- unlist(test$p.value)
    t.df$alternative <- unlist(test$alternative)
    t.df$method <- unlist(test$method)
    t.df <- as.data.frame(t.df) %>% filter(p.value < 0.05,alternative=="greater")
    write.csv(t.df,paste0(positive[1],"_vs_",negative[n],".csv"))   #输出csv文件
    y <- t.df$microbe
    if(n==1){
      x <- t.df$microbe
    } else {
      x <- intersect(x,y)
      print(length(x))
    }
    
    #sink(paste0(positive[1],"_vs_",negative[n],".csv"))
    #as.data.frame(test) %>% filter(p.value < 0.05)
    #sink()
    n=n+1
  }
  #x即为positive相对于所有negative具有秩和检验显著性的特征
  print(x)
  write.csv(x,paste0(positive[1],"_specific_feature.csv"))
  
  #plot 需要手动修改
  candidate <- x
  counts <- read.csv("counts-G_norm_decontamed.txt",sep="\t",header = T, row.names = 1)
  {
    j=1
    {
      test={}
      while(j<=length(candidate)){
        target <- candidate[j]
        type <- c("CRC","NC_PKU","STAD","HCC","LUAD","ESCA")
        i=1
        gene={}
        while(i<=length(type)){
          temp <- counts[grep(target,rownames(counts),fixed=TRUE)[1],grep(type[i],colnames(counts),fixed=TRUE)]
          temp.t <- t(temp)
          colnames(temp.t) <- c("counts")
          temp.df <- as.data.frame(temp.t)
          #k=1
          #OutVals <- boxplot(temp.t,plot=FALSE)$out
          #OutVals <- unique(OutVals)
          #OutVals <- sort(OutVals,decreasing = T)
          #while(k<=length(OutVals)){
          #  # return integer(0), temp.df <- temp.df
          #  if(length(grep(OutVals[k],temp.df$counts))==0){
          #    temp.df <- temp.df
          #  } else {
          #    temp.df <- as.data.frame(temp.df[-grep(OutVals[k],temp.df$counts),])
          #    colnames(temp.df) <- "counts"
          #    #print(length(temp.df$counts))
          #  }
          #  k=k+1
          #}
          temp.df$CancerType<-rep(type[i],times=length(temp.df$counts))
          temp.df$microbe<-rep(target,times=length(temp.df$counts))
          gene <- rbind(gene,temp.df)
          i=i+1
        }
        
        
        my_comparisons <- list(c("CRC","NC_PKU"),c("CRC", "STAD"),c("CRC","HCC"), c("CRC", "LUAD"), c("CRC", "ESCA"))
        gene$CancerType <- factor(gene$CancerType,levels=type)
        p <- ggplot(gene,aes(x=CancerType,y=counts,fill=CancerType))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
          geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
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
            axis.title.y = element_text(face="bold",color="black", size=24))
        
        m=10
        if(m>0){
          p <- p+stat_compare_means(comparisons = my_comparisons,
                                    method = "wilcox.test",
                                    method.args = list(alternative = "greater"),
                                    label = "p.signif"
          )+labs(x="",y="Relative Abundance (%)",title=paste0(target,"\nwilcox.test.greater"), face="bold",fill="Batch",width=10,height=8)
          ggsave(p,filename = paste0(target,"_greater.pdf"),path = "./STAD_filtered_greater/remove_outlier")
        } else {
          p <- p+stat_compare_means(comparisons = my_comparisons,
                                    method = "wilcox.test",
                                    method.args = list(alternative = "less"),
                                    label = "p.signif"
          )+labs(x="",y="Relative Abundance (%)",title=paste0(target,"\nwilcox.test.less"), face="bold",fill="Batch",width=10,height=8)
          ggsave(p,filename = paste0(target,"_less.pdf"),path="./")
        }
        j=j+1
      }
    }
  }
  
}

#### ribosomal protein mRNA selection workflow
library(clusterProfiler) ## enrichment
library(biomaRt) ## ID convert
library(org.Hs.eg.db) ## annotation for human
library(enrichplot) ##gseaplot
#pathway gene accession
{
hsa_kegg <- clusterProfiler::download_KEGG("hsa")
names(hsa_kegg)
PATH2NAME <- hsa_kegg$KEGGPATHID2NAME
PATH2ID <- hsa_kegg$KEGGPATHID2EXTID
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by="from")
colnames(PATH_ID_NAME) <- c("KEGGID", "ENTREZID", "DESCRPTION")
#ID conversion
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
entrezgene <- PATH_ID_NAME$ENTREZID
ensembl_gene_id<- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                        filters = "entrezgene_id",
                        values=entrezgene, mart= mart,useCache = FALSE)
PATH_ID_NAME <- merge(PATH_ID_NAME, ensembl_gene_id, by.x= "ENTREZID",by.y= "entrezgene_id")
Ribosome <- subset(PATH_ID_NAME,DESCRPTION=="Ribosome")
Ribosome_biogensis <- subset(PATH_ID_NAME,DESCRPTION=="Ribosome biogenesis in eukaryotes")
Sulfur_metabolism <- subset(PATH_ID_NAME,DESCRPTION=="Sulfur metabolism")
Malaria <- subset(PATH_ID_NAME,DESCRPTION=="Malaria")
kappa <- subset(PATH_ID_NAME,DESCRPTION=="NF-kappa B signaling pathway")
IL17 <- subset(PATH_ID_NAME,DESCRPTION=="IL-17 signaling pathway")
}

## boxpot for one pathway(need to manually specify "less","greater" or "two.sided")
#all comparison take healthy control as positive, which means "greater" is greater in healthy, "less" is less in healthy. 
{
  test_direction <- "less"
  pathway_name <- "Ribosome"
  pathway <- Ribosome
{ 
  setwd(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/",pathway_name,"/boxplot/THCA"))
  #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/GSE68086_TPM.txt",sep="\t",header = T, row.names = 1)
  #positive <- c("HD")
  #negative <- c("Breast","Liver","CRC","Lung","Panc","NSCLC","GBM")
  #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/PBMC_TPM.txt",sep="\t",header = T, row.names = 1)
  #positive <- c("NC")
  #negative <- c("CRC")
  counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/WCHSU_PTC_16_TPM.txt",sep="\t",header = T, row.names = 1)
  positive <- c("N")
  negative <- c("THCA")
  
  
  candidate <- pathway$ensembl_gene_id
  candidate <- unique(candidate)
  n=1
  counts = log2(counts+1)
  while(n<=length(negative)){
    #logFC <- apply(matrix_log2counts[,grep(positive[1],colnames(matrix_log2counts),fixed=TRUE)], 1, mean)-apply(matrix_log2counts[,grep(negative[n],colnames(matrix_log2counts),fixed=TRUE)], 1, mean)
    #logFC <- as.data.frame(logFC)
    j=1
    #while(j<=length(candidate)){
    #  print(grep(candidate[j],rownames(counts),fixed=TRUE))
    #  j=j+1
    #}
    test={}
    type <- c(positive[1],negative[n])
    while(j<=length(candidate)){
      target <- candidate[j]
      if(length(grep(target,rownames(counts)))==0) {
        #print(paste0(target," is rRNA?"))
        j=j+1
      } else {
      i=1
      gene={}
      while(i<=length(type)){
        temp <- counts[grep(target,rownames(counts),fixed=TRUE)[1],grep(type[i],colnames(counts),fixed=TRUE)]
        temp.t <- t(temp)
        colnames(temp.t) <- c("counts")
        temp.df <- as.data.frame(temp.t)
        #remove outlier
        #k=1
        #OutVals <- boxplot(temp.t,plot=FALSE)$out
        #OutVals <- unique(OutVals)
        #OutVals <- sort(OutVals,decreasing = T)
        #while(k<=length(OutVals)){
        #  # return integer(0), temp.df <- temp.df
        #  if(length(grep(OutVals[k],temp.df$counts))==0){
        #    temp.df <- temp.df
        #  } else if(OutVals[k]==0) {
        #    temp.df <- temp.df
        #  }
        #  else {
        #    temp.df <- as.data.frame(temp.df[-grep(OutVals[k],temp.df$counts,fixed=TRUE),])
        #    colnames(temp.df) <- "counts"
        #    #print(length(temp.df$counts))
        #  }
        #  k=k+1
        #}
        temp.df$CancerType<-rep(type[i],times=length(temp.df$counts))
        temp.df$gene<-rep(target,times=length(temp.df$counts))
        gene <- rbind(gene,temp.df)
        i=i+1
      }
      trend <- mean(subset(gene,CancerType==positive[1])$counts)-mean(subset(gene,CancerType==negative[n])$counts)
      #plot
      {
      #  my_comparisons <- type 
      #  gene$CancerType <- factor(gene$CancerType,levels=type)
      #  p <- ggplot(gene,aes(x=CancerType,y=counts,fill=CancerType))+ 
      #  geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      #  geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+ 
      #  scale_fill_brewer(palette="Blues") + 
      #  theme_bw()+
      #  theme(legend.position="right",
      #    panel.grid=element_blank(),
      #    panel.border=element_blank(),
      #    axis.line = element_line(size=1, colour = "black"),
      #    legend.title = element_text(face="bold", color="black",family = "Arial", size=10),
      #    legend.text= element_text(face="bold", color="black",family = "Arial", size=10),
      #    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      #    axis.text.x = element_text(face="bold", color="black", size=10,angle = 45,hjust = 1),
      #    axis.text.y = element_text(face="bold",  color="black", size=10),
      #    axis.title.x = element_text(face="bold", color="black", size=10),
      #    axis.title.y = element_text(face="bold",color="black", size=10))
      }
      #wilcox.test       
      if(trend>0){
        #p <- p+stat_compare_means(comparisons = my_comparisons,
        #                          method = "wilcox.test",
        #                          method.args = list(alternative = "greater"),
        #                          label = "p.signif")+labs(x="",y="relative abundance",title=paste0(target,"_wilcox.test.greater"), face="bold")
        #ggsave(p,filename = paste0(target,"_greater.pdf"),path = "./")
        t <- wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="greater" )
      } else if(trend<0){
        #p <- p+stat_compare_means(comparisons = my_comparisons,
        #                          method = "wilcox.test",
        #                          method.args = list(alternative = "less"),
        #                          label = "p.signif")+labs(x="",y="relative abundance",title=paste0(target,"_wilcox.test.less"), face="bold")
        #ggsave(p,filename = paste0(target,"_less.pdf"),path="./")
        t <- wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="less" )
      } else {
        t <- wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="two.sided" )
      }
      t$gene <- rep(target,times=length(t$p.value)) 
      t$positive_vs_negative <- rep(paste0(positive[1],"_vs_",negative[n]),times=length(t$p.value))
      test <- rbind(test,t)
      j=j+1
      }
    }
    t.df={}
    test <- as.data.frame(test)
    t.df$gene <- unlist(test$gene)
    t.df$positive_vs_negative <- unlist(test$positive_vs_negative)
    t.df$p.value <- unlist(test$p.value)
    t.df$alternative <- unlist(test$alternative)
    t.df$method <- unlist(test$method)
    t.df <- as.data.frame(t.df) %>% filter(p.value < 0.05,alternative==test_direction)
    write.csv(t.df,paste0(positive[1],"_vs_",negative[n],".csv"))   #输出csv文件
    y <- t.df$gene
    if(n==1){
      x <- t.df$gene
    } else {
      x <- intersect(x,y)
      print(length(x))
    }
    
    #sink(paste0(positive[1],"_vs_",negative[n],".csv"))
    #as.data.frame(test) %>% filter(p.value < 0.05)
    #sink()
    n=n+1
  }
  #x即为positive相对于所有negative具有秩和检验显著性的特征
  print(x)
  write.csv(x,paste0(positive[1],"_specific_feature.csv"))
}

#boxplot
#x <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/Ribosome.csv") 
#x <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/N_vs_PDAC.csv") 
candidate <- x
#candidate <-c("ENSG00000148303","RPL7A","ENSG00000169100","SLC25A6","NAPIL1","RACK1","DEK","RPS10","APPR19","HMGB1","ADAM","EEF","EIF","RPS10-NUD")
#candidate <- c("ENSG00000167526","ENSG00000174748")
#candidate <- c("ENSG00000167526","ENSG00000174748","ENSG00000136942")
{
  j=1
  {
    test={}
    while(j<=length(as.matrix(candidate))){
      target <- as.matrix(candidate)[j]
      #type <- c("NC_PKU","CRC","STAD","LUAD","ESCA_PKU") #pico_PKU
      #type <- c("NC_ShH","HCC_ShH") #pico_ShH
      #type <- c("NC_ChQ","HCC_ChQ") #pico_ChQ
      #type <- c("NC","CRC","CHD","HCC","PAAD") #exorbase
      #type <- c("N","PDAC","CP") #gut GSE133684
      #type <- c("HD","NSCLC") # GSE89843 all batch
      #type <- c("Vumc.HD","Vumc.NSCLC") # GSE89843 in vumc
      #type <- c("HD","Breast","Liver","CRC","Lung","Panc","NSCLC","GBM") # GSE68086
      #type <- c("N","T","plasma","HD","PBMC") #pico_PBMC
      #type <- c("N","HCC") # GSE120663
      #type <- c("NC","CRC") #PBMC
      type <- c("N","THCA") #THCA
      i=1
      gene={}
      while(i<=length(type)){
        temp <- counts[grep(target,rownames(counts),fixed=TRUE)[1],grep(type[i],colnames(counts),fixed=TRUE)]
        temp.t <- t(temp)
        colnames(temp.t) <- c("counts")
        temp.df <- as.data.frame(temp.t)
        #k=1
        #OutVals <- boxplot(temp.t,plot=FALSE)$out
        #OutVals <- unique(OutVals)
        #OutVals <- sort(OutVals,decreasing = T)
        #while(k<=length(OutVals)){
        #  # return integer(0), temp.df <- temp.df
        #  if(length(grep(OutVals[k],temp.df$counts))==0){
        #    temp.df <- temp.df
        #  } else {
        #    temp.df <- as.data.frame(temp.df[-grep(OutVals[k],temp.df$counts),])
        #    colnames(temp.df) <- "counts"
        #    #print(length(temp.df$counts))
        #  }
        #  k=k+1
        #}
        temp.df$CancerType<-rep(type[i],times=length(temp.df$counts))
        temp.df$microbe<-rep(target,times=length(temp.df$counts))
        gene <- rbind(gene,temp.df)
        i=i+1
      }
      
      
      #my_comparisons <- list(c("NC_PKU","CRC"),c("NC_PKU","STAD"),c("NC_PKU","LUAD"),c("NC_PKU","ESCA_PKU")) #pico_PKU
      #my_comparisons <- list(c("NC_ShH","HCC_ShH")) #pico_ShH
      #my_comparisons <- list(c("NC_ChQ","HCC_ChQ")) #pico_ChQ
      #my_comparisons <- list(c("NC","CRC"),c("NC","CHD"),c("NC","HCC"),c("NC","PAAD")) #exorbase
      #my_comparisons <- list(c("N","PDAC"),c("N","CP")) #gut
      #my_comparisons <- list(c("HD","NSCLC")) #GSE89843 platelet all batch
      #my_comparisons <- list(c("Vumc.HD","Vumc.NSCLC")) #GSE89843 in vumc
      #my_comparisons <- list(c("HD","Breast"),c("HD","Liver"),c("HD","CRC"),c("HD","Lung"),c("HD","Panc"),c("HD","NSCLC"),c("HD","GBM"))  #GSE68086
      #my_comparisons <- list(c("N","T"),c("plasma","HD"),c("plasma","PBMC")) #pico_PBMC
      #my_comparisons <- list(c("N","HCC")) # GSE120663
      #my_comparisons <- list(c("NC","CRC")) # PBMC
      my_comparisons <- list(c("N","THCA")) # THCA
      gene$CancerType <- factor(gene$CancerType,levels=type)
      p <- ggplot(gene,aes(x=CancerType,y=counts,fill=CancerType))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
        geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
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
          axis.title.y = element_text(face="bold",color="black", size=24))
      
      m=test_direction
      if(m=="greater"){
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  method.args = list(alternative = "greater"),
                                  label = "p.signif"
        )+labs(x="",y="TPM",title=paste0(target,"\nwilcox.test.greater"), face="bold",fill="Type",width=10,height=8)
        ggsave(p,filename = paste0(target,"_greater.pdf"),path = "./")
      } else if(m=="less") {
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  method.args = list(alternative = "less"),
                                  label = "p.signif"
        )+labs(x="",y="TPM",title=paste0(target,"\nwilcox.test.less"), face="bold",fill="Type",width=10,height=8)
        ggsave(p,filename = paste0(target,"_less.pdf"),path="./")
      } else if(m=="two.sided") {
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  method.args = list(alternative = "two.sided"),
                                  label = "p.signif"
        )+labs(x="",y="TPM",title=paste0(target,"\nwilcox.test.two.sided"), face="bold",fill="Type",width=10,height=8)
        ggsave(p,filename = paste0(target,"_less.pdf"),path="./")
      } else {
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  method.args = list(alternative = "two.sided"),
                                  label = "p.signif"
        )+labs(x="",y="TPM",title=paste0(target,"\nwilcox.test.two.sided"), face="bold",fill="Type",width=10,height=8)
        ggsave(p,filename = paste0(target,"_less.pdf"),path="./")
      }
      j=j+1
    }
  }
}

}

# heatmap
{
target <- unique(Ribosome$ensembl_gene_id)
#pico
#counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/TPM_all_pico.txt",sep="\t",header = T, row.names = 1)
#type <- c("NC_ShH","HCC_ShH")
#type <- c("NC_ChQ","HCC_ChQ")
#type <- c("NC_PKU","CRC","STAD","LUAD","ESCA_PKU")
#exoRbase
#type <- c("NC","CRC","CHD","HCC","PAAD")
#counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/exoRbase_TPM.txt",sep="\t",header = T, row.names = 1)
#GSE133684
#type <- c("N","CP","PDAC")
#counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/GSE133684_exp_TPM-all.txt",sep="\t",header = T, row.names = 1)
#GSE89843 vumc
#type <- c("Vumc.HD","Vumc.NSCLC")
#counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/GSE89843_TPM.txt",sep="\t",header = T, row.names = 1)
#GSE68086
#type <- c("HD","Breast","Liver","CRC","Lung","Panc","NSCLC","GBM")
#counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/GSE68086_TPM.txt",sep="\t",header = T, row.names = 1)
#PBMC
#type <- c("NC","CRC")
#counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/PBMC_TPM.txt",sep="\t",header = T, row.names = 1)
#THCA
type <- c("N","THCA")
counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/WCHSU_PTC_16_TPM.txt",header = T, row.names = 1,sep="\t")

df.type=as.data.frame(array(dim=c(1,1))) #设置数据框行列
i=1
while(i<=length(type)){
  df <- counts[,grep(type[i],colnames(counts),fixed=TRUE)]
  df.type <- cbind(df.type,df)
  i=i+1
}
df.type <- df.type[,-1]

# progress
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = length(target), clear = FALSE, width= 60) 
df.target={}
j=1
while(j<=length(target)){
  df <- df.type[grep(target[j],rownames(df.type),fixed=TRUE),]
  df.target <- rbind(df.target,df)
  j=j+1
  pb$tick()
  Sys.sleep(1 / 100)
}

heatmap <- df.target
#colnames(heatmap) <- gsub("HDControl","HD",colnames(heatmap)) # for GSE89843 vumc

#for inhouse processed datasets, necessary to extract ENSEMBL ID, for downloaded processed dataset, not necessary
rownames(heatmap) <- as.character(lapply(strsplit(rownames(df.target),"\\|"),function(x) x[1])) #1 for ENSG at each line begin, 3 for pico(ENSG at the third col)
rownames(heatmap) <- as.character(lapply(strsplit(rownames(heatmap),"\\."),function(x) x[1])) #1 for ENSG at each line begin, 3 for pico(ENSG at the third col)

# key type transvert
symbol <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(heatmap), columns=c("SYMBOL"), keytype="ENSEMBL")
symbol$ENSEMBL #check and find 77st row ensembl ID repeatly show up(for gut GSE133864) #rownames(heatmap) <- symbol[-78,]$SYMBOL
#rownames(heatmap) <- symbol[-33,]$SYMBOL #check and find 32st row ensembl ID repeatly show up(for GSE89843)
#rownames(heatmap) <- symbol[-76,]$SYMBOL #check and find 76st row ensembl ID repeatly show up(for GSE68086)
#rownames(heatmap) <- symbol[-119,]$SYMBOL #check and find 119st row ensembl ID repeatly show up(for pico)
rownames(heatmap) <- symbol[-136,]$SYMBOL #check and find 119st row ensembl ID repeatly show up(for THCA)

column_annotation = as.data.frame(array(dim=c(length(colnames(heatmap)),1))) 
rownames(column_annotation) <- colnames(heatmap)
#GSE68086 must manully revise
#write.csv(colnames(heatmap),"heatmap_annotation.csv")
#column_annotation <- read.csv("heatmap_annotation.csv",row.names = 1,header = F)
column_annotation$V1 <- as.character(lapply(strsplit(colnames(heatmap),".",fixed = TRUE),function(x) x[2])) ##exoRBase "_",x[1]; pico ".",x[1];GSE89843, ".", x[2]
colnames(column_annotation) <- c("Type")



pheatmap(
  heatmap,
  annotation_col = column_annotation,
  scale = "row",
  cluster_cols = FALSE,cluster_rows = FALSE,
  show_colnames=FALSE, 
  #cluster_cols = FALSE,
  colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 5)
}

# GSEA (included gene rank generation by edger_glmlrt)
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome")
  mat_raw <- read.table("Raw_counts_passed_QC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
  #mat_raw <- read.table("exoRbase_rawcounts.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') 
  #mat_raw <- read.table("GSE68086_TEP_data_matrix_forDE.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
  #mat_raw <- read.table("GSE89843_TEP_Count_Matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
  #mat_raw <- read.table("PBMC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
  #des <- read.csv("des_HCC_ShHvsNC_ShH.csv", header = TRUE, check.names=FALSE, sep=',')  #读入数据
  des <- read.csv("des_CancervsNC_PKU.csv", header = TRUE, check.names=FALSE, sep=',')
  #des <- read.csv("des_HCC_ChQvsNC_ChQ.csv", header = TRUE, check.names=FALSE, sep=',')
  #des <- read.csv("des_exoRbase.csv", header = TRUE, check.names=FALSE, sep=',')
  #des <- read.csv("des_GSE68086.csv", header = TRUE, check.names=FALSE, sep=',')
  #des <- read.csv("des_GSE89843.csv", header = TRUE, check.names=FALSE, sep=',')
  #des <- read.csv("des_PBMC.csv", header = TRUE, check.names=FALSE, sep=',')
  
  samples <- des$samples
  #samples <- gsub(".","-",samples,fixed=TRUE)
  #samples <- gsub("X","",samples,fixed=TRUE)
  group <- des$group
  batch <- des$batch
  
  i=1
  mat=as.data.frame(array(dim=c(length(rownames(mat_raw)),1)))
  while (i<=length(samples)) {
    temp <- mat_raw[,which(colnames(mat_raw)==samples[i])]
    #temp <- mat_raw[,grep(samples[i],colnames(mat_raw),fixed=TRUE)]
    #print(i)
    #print(which(colnames(mat_raw)==samples[i]))
    temp <- as.data.frame(temp)
    colnames(temp) <- samples[i]
    mat <- cbind(mat,temp)
    i=i+1
  }
  
  mat <- mat[,-1]
  rownames(mat) <- rownames(mat_raw)
  
  y <- DGEList(counts=mat, samples=samples, group=group)     #规定输入格式，必须包含这几个参数
  y <- calcNormFactors(y, method="TMM")                      #归一化处理
  y <- estimateDisp(y)
  test <- exactTest(y, pair = c("negative","positive"), dispersion = "auto")
  #广义线性模型可以针对多因素的情况，在这个场景下，只比较cancer与normal用exactTest就足够了
  {
    #design <- model.matrix(~0+group)   
    #y <- estimateDisp(y, design)          #广义线性模型计算离散度（common&trended&tagwise）
    #fit <- glmFit(y, design)              #差异表达分析函数                     
    #test <- glmLRT(fit, coef=2)           #差异表达分析函数 
  }
  
  
  
  res <- topTags(test, n=nrow(mat), sort.by='none')     #输出计算差异基因的结果
  res <- cbind(res$table, baseMean=2^(res$table$logCPM)) #输出差异基因结果加上baseMean一栏
  mapped_names <- colnames(res)                   #后面都是整理+输出
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  write.table(res, "edger_exact.csv", sep=',', quote=FALSE, row.names=TRUE) #输出文件，规定名字
}

#make ranked genelist
{
res.sort <- res[sort(res$log2FoldChange,decreasing = T, index.return = T)$ix,]
res.sort <- res.sort[grep("ENSG",rownames(res.sort)),]

genelist <- as.numeric(res.sort$log2FoldChange)
names(genelist) <- as.character(lapply(strsplit(rownames(res.sort),"\\|"), function(x) x[1]))
names(genelist) <- as.character(lapply(strsplit(names(genelist),".",fixed = TRUE), function(x) x[1]))

}

#make geneset from KEGG
{
  hsa_kegg <- clusterProfiler::download_KEGG("hsa")
  names(hsa_kegg)
  PATH2NAME <- hsa_kegg$KEGGPATHID2NAME
  PATH2ID <- hsa_kegg$KEGGPATHID2EXTID
  PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by="from")
  colnames(PATH_ID_NAME) <- c("KEGGID", "ENTREZID", "DESCRPTION")
  #ID conversion
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  entrezgene <- PATH_ID_NAME$ENTREZID
  ensembl_gene_id<- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                          filters = "entrezgene_id",
                          values=entrezgene, mart= mart,useCache = FALSE)
  PATH_ID_NAME <- merge(PATH_ID_NAME, ensembl_gene_id, by.x= "ENTREZID",by.y= "entrezgene_id")
  ensg <- PATH_ID_NAME$ensembl_gene_id
  gene_symbol <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),
                       filters = "ensembl_gene_id",
                       values=ensg, mart= mart,useCache = FALSE)
  PATH_ID_NAME <- merge(PATH_ID_NAME, gene_symbol, by.x= "ensembl_gene_id",by.y= "ensembl_gene_id")
  PATH_ID_NAME <- unique(PATH_ID_NAME)
  #gmt <- PATH_ID_NAME[,3:4] #from web
  gmt <- PATH_ID_NAME[,c(4,1)] # local
  colnames(gmt) <- c("ont","gene")
}

#GSEA by clusterprofile and plot
{
egmt2 <- GSEA(genelist,minGSSize = 0, maxGSSize = 1000,
              pvalueCutoff = 1, TERM2GENE = gmt[which(gmt$ont=="Ribosome"),])

gsea <- gseaplot2(egmt2,1)
pdf("Ribosome_PBMC_CRCvsNC.pdf")
plot(gsea)
dev.off()
}

#if cannot get enrichment results, use local GSEA to plot
rnk <- as.data.frame(as.matrix(genelist))
colnames(rnk) <- c("#generank")
write.table(rnk,"./GSEA/rank_PBMC.csv",sep = "\t")
# need to manually split by tab

#pathway level boxplot (like gepia provided)

#summary pathway level matrix
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/")
#counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/GSE136243_gene_TPM_all_samples.txt",sep="\t",header = T, row.names = 1)
#counts <- read.csv("/Users/yuhuan/Desktop/PBMC/GSE120663_TPM.csv",header = TRUE, row.names = 1)
#counts <- counts[,-1] # for GSE40174, whose 1st colum is ENTREZGENE_ID
#counts <- counts[,-c(1,2)] # for GSE136243, whose 1st colum is gene_name, 2nd colum is gene_biotype
#counts <- read.csv("/Users/yuhuan/Desktop/PBMC/pico_PBMC_TPM.csv",header = TRUE, row.names = 1)
#counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/PBMC_TPM.txt",sep="\t",header = T, row.names = 1)
#counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/WCHSU_PTC_16_TPM.txt",sep="\t",header = T, row.names = 1)
#counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/pico_tissue_TPM.txt",sep="\t",header = T, row.names = 1)

#20210611 multiomics plasma
counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/plasma/TPM_intron_spanning_noMTRNA.txt",sep="\t",header = T, row.names = 1)
counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/tissue/tissue_TPM_intron_spanning_noMTRNA.txt",sep="\t",header = T, row.names = 1)
counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/05.Deconvolution/PBMC/TPM_intron_spanning_noMTRNA.txt",sep="\t",header = T, row.names = 1)

counts <- log2(as.matrix(counts)+1)

# set progress
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = length(unique(PATH_ID_NAME$DESCRPTION)), clear = FALSE, width= 60) 
#get all genes in each pathway in KEGG
n=1
pathway_count={}
while(n <= length(unique(PATH_ID_NAME$DESCRPTION))){
pathway_name <- as.character(unique(PATH_ID_NAME$DESCRPTION)[n])
pathway <- subset(PATH_ID_NAME,DESCRPTION==pathway_name)

candidate <- pathway$ensembl_gene_id #for colnames with ensemble id
#candidate <- pathway$hgnc_symbol #for colnames with gene symbol

candidate <- unique(candidate)

#summary one pathway total count in each sample and make 1 row dataframe
j=1
pathway_gene_count={}
while(j<=length(candidate)){
  target <- candidate[j]
  if(length(grep(target,rownames(counts)))==0) {
    #print(paste0("No ",target," in this dataset."))
    j=j+1
  } else {
    #temp <- counts[which(rownames(counts)==target),]
    temp <- counts[grep(target,rownames(counts),fixed=TRUE),]
    pathway_gene_count <- rbind(pathway_gene_count,temp)
    j=j+1
  }
}
if(is.null(pathway_gene_count)){
  n=n+1
  next;
  } else {
one_pathway <- as.data.frame(t(colMeans(pathway_gene_count))) #colSums
rownames(one_pathway) <- pathway_name
pathway_count <- rbind(pathway_count,one_pathway)
n=n+1
}
pb$tick()
Sys.sleep(1 / 100)

}
write.csv(pathway_count,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/PBMC_TPM_noMT_pathway_level_log2_mean.csv")
}

#pathway_count <- read.csv("./GSE89843/GSE89843_TPM_pathway_level.csv",header = TRUE , row.names = 1)
#pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/PBMC/PBMC_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/pico/pico_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/pico/pico_TPM_pathway_level_log2_gender.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/GSE27562_RNA_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/pico_tissue/pico_tissue_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/PBMC/PBMC_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE68086/GSE68086_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE89843/GSE89843_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/exoRBase/exoRBase_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE133684/GSE133684_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE120663/GSE120663_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE40174/GSE40174_normalized_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/THCA/THCA_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/tissue_TPM_noMT_pathway_level_log2_mean.csv",header = TRUE , row.names = 1)
pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/plasma_TPM_noMT_pathway_level_log2_mean.csv",header = TRUE , row.names = 1)


# filter significant changed pathway by wilcox test
{ 
  #GSE68086
  #positive <- c("HD")
  #negative <- c("Breast","Liver","CRC","Lung","Panc","NSCLC","GBM")
  #pico_PKU
  #positive <- c("NC_PKU")
  #negative <- c("CRC","STAD","LUAD","ESCA_PKU")
  #pico_ShH
  #positive <- c("NC_ShH")
  #negative <- c("HCC_ShH")
  #exoRBase
  #positive <- c("NC")
  #negative <- c("CRC","CHD","HCC","PAAD")
  #GSE89843
  #positive <- c("Vumc.HD")
  #negative <- c("Vumc.NSCLC")
  #GSE133684
  #positive <- c("N")
  #negative <- c("PDAC")
  #GSE40174
  #positive <- c("HD")
  #negative <- c("PDAC")
  #GSE131512
  #positive <- c("N")
  #negative <- c("C")
  #GSE136243
  #positive <- c("N")
  #negative <- c("AD")
  #PBMC
  #positive <- c("NC")
  #negative <- c("CRC")
  #THCA
  #positive <- c("N")
  #negative <- c("THCA")
  #stage_pico
  #positive <- c("NC_PKU")
  #negative <- c("I.CRC","II.CRC","III.CRC","IV.CRC","I.STAD","II.STAD","III.STAD","IV.STAD")
  #GSE27562
  #positive <- c("normal")
  #negative <- c("benign","malignant","post.surgery","ectopic")
  #multi-omics_tissue
  #positive <- c("N")
  #negative <- c("T")
  #multi-omics_plasma
  positive <- c("NC")
  negative <- c("CRC")
  
  candidate <- rownames(pathway_count)
  candidate <- unique(candidate)
  candidate <- "Ribosome"
  n=1
  while(n<=length(negative)){
    j=1
    test={}
    type <- c(positive[1],negative[n])
    while(j<=length(candidate)){
      target <- candidate[j]
      if(length(grep(target,rownames(pathway_count)))==0) {
        #print(paste0("No ",target," in this dataset."))
        j=j+1
      } else {
        i=1
        gene={}
        while(i<=length(type)){
          temp <- pathway_count[which(rownames(pathway_count)==target),grep(type[i],colnames(pathway_count),fixed=TRUE)]
          temp.t <- t(temp)
          colnames(temp.t) <- c("counts")
          temp.df <- as.data.frame(temp.t)
          temp.df$CancerType<-rep(type[i],times=length(temp.df$counts))
          temp.df$gene<-rep(target,times=length(temp.df$counts))
          gene <- rbind(gene,temp.df)
          i=i+1
        }
        trend <- mean(subset(gene,CancerType==positive[1])$counts)-mean(subset(gene,CancerType==negative[n])$counts)
        #plot
        {
          #  my_comparisons <- type 
          #  gene$CancerType <- factor(gene$CancerType,levels=type)
          #  p <- ggplot(gene,aes(x=CancerType,y=counts,fill=CancerType))+ 
          #  geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
          #  geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+ 
          #  scale_fill_brewer(palette="Blues") + 
          #  theme_bw()+
          #  theme(legend.position="right",
          #    panel.grid=element_blank(),
          #    panel.border=element_blank(),
          #    axis.line = element_line(size=1, colour = "black"),
          #    legend.title = element_text(face="bold", color="black",family = "Arial", size=10),
          #    legend.text= element_text(face="bold", color="black",family = "Arial", size=10),
          #    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
          #    axis.text.x = element_text(face="bold", color="black", size=10,angle = 45,hjust = 1),
          #    axis.text.y = element_text(face="bold",  color="black", size=10),
          #    axis.title.x = element_text(face="bold", color="black", size=10),
          #    axis.title.y = element_text(face="bold",color="black", size=10))
        }
        #wilcox.test       
        if(trend>0){
          #p <- p+stat_compare_means(comparisons = my_comparisons,
          #                          method = "wilcox.test",
          #                          method.args = list(alternative = "greater"),
          #                          label = "p.signif")+labs(x="",y="relative abundance",title=paste0(target,"_wilcox.test.greater"), face="bold")
          #ggsave(p,filename = paste0(target,"_greater.pdf"),path = "./")
          t <- wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="greater" )
        } else if(trend<0){
          #p <- p+stat_compare_means(comparisons = my_comparisons,
          #                          method = "wilcox.test",
          #                          method.args = list(alternative = "less"),
          #                          label = "p.signif")+labs(x="",y="relative abundance",title=paste0(target,"_wilcox.test.less"), face="bold")
          #ggsave(p,filename = paste0(target,"_less.pdf"),path="./")
          t <- wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="less" )
        } else {
          t <- wilcox.test(subset(gene,CancerType==positive[1])$counts,subset(gene,CancerType==negative[n])$counts,alternative="two.sided" )
        }
        t$gene <- rep(target,times=length(t$p.value)) 
        t$positive_vs_negative <- rep(paste0(positive[1],"_vs_",negative[n]),times=length(t$p.value))
        test <- rbind(test,t)
        j=j+1
      }
    }
    t.df={}
    test <- as.data.frame(test)
    t.df$gene <- unlist(test$gene)
    t.df$positive_vs_negative <- unlist(test$positive_vs_negative)
    t.df$p.value <- unlist(test$p.value)
    t.df$alternative <- unlist(test$alternative)
    t.df$method <- unlist(test$method)
    t.df <- as.data.frame(t.df) %>% filter(p.value < 0.05,alternative=="less") # 在filter()中添加alternative=="greater"可以特异的挑出在正常人中上调的通路,alternative=="less"可以特异的挑出在正常人中下调的通路
    write.csv(t.df,paste0(positive[1],"_vs_",negative[n],".pathway_level.csv"))   #输出csv文件
    y <- t.df$gene
    if(n==1){
      x <- t.df$gene
    } else {
      x <- intersect(x,y)
      print(length(x))
    }
    
    #sink(paste0(positive[1],"_vs_",negative[n],".csv"))
    #as.data.frame(test) %>% filter(p.value < 0.05)
    #sink()
    n=n+1
  }
  #x即为positive相对于所有negative具有秩和检验显著性的特征
  print(x)
  write.csv(x,paste0(positive[1],"_specific_feature.pathway_level.csv"))
}

#boxplot
#x <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/GSE89843/Upregulated in Cancer/Vumc.HD_specific_feature.pathway_level.csv",header = FALSE) 
candidate <- x

#candidate <-c("ENSG00000148303","RPL7A","ENSG00000169100","SLC25A6","NAPIL1","RACK1","DEK","RPS10","APPR19","HMGB1","ADAM","EEF","EIF","RPS10-NUD")
#candidate <- c("ENSG00000167526","ENSG00000174748")
#candidate <- c("ENSG00000167526","ENSG00000174748","ENSG00000136942")
{
  j=1
  {
    test={}
    while(j<=length(as.matrix(candidate))){
      target <- as.matrix(candidate)[j]
      #type <- c("NC_PKU","CRC","STAD","LUAD","ESCA_PKU","NC_ShH","HCC_ShH","NC_ChQ","HCC_ChQ") #pico
      #type <- c("NC_PKU","CRC","STAD","LUAD","ESCA_PKU")#pico_PKU
      #type <- c("NC_ShH","HCC_ShH")#pico_ShH
      #type <- c("NC","CRC","CHD","HCC","PAAD") #exorbase
      #type <- c("N","PDAC") #gut GSE133684
      #type <- c("HD","NSCLC") # GSE89843 all batch
      #type <- c("Vumc.HD","Vumc.NSCLC") # GSE89843 in vumc
      #type <- c("HD","Breast","Liver","CRC","Lung","Panc","NSCLC","GBM") # GSE68086
      #type <- c("HD","PDAC") #GSE40174
      #type <- c("N","C") #GSE131512
      #type <- c("N","AD") #GSE136243 
      #type <- c("N","HCC") #GSE120663
      #type <- c("NC","CRC") #PBMC
      #type <- c("N","THCA") #THCA
      #type <- c("NC_PKU","I.CRC","II.CRC","III.CRC","IV.CRC","I.STAD","II.STAD","III.STAD","IV.STAD","I.LUAD","II.LUAD","III.LUAD","IV.LUAD",
      #          "I.ESCA_PKU","II.ESCA_PKU","III.ESCA_PKU","IV.ESCA_PKU","NC_ShH","I.HCC_ShH","II.HCC_ShH","III.HCC_ShH","IV.HCC_ShH") #stage
      #type <- c("M.NC_PKU","F.NC_PKU","M.CRC","F.CRC","M.STAD","F.STAD","M.LUAD","F.LUAD","M.ESCA_PKU","F.ESCA_PKU","M.HCC_ShH","F.HCC_ShH") #gender
      #type <- c("normal","benign","malignant","post.surgery","ectopic")
      #type <- c("N","T") #pico_tissue
      #type <- c("NC","CRC","STAD")
      #type <- c("N","T")
      type <- c("NC","CRC")
      
      i=1
      gene={}
      while(i<=length(type)){
        temp <- pathway_count[which(rownames(pathway_count)==target),grep(type[i],colnames(pathway_count),fixed=TRUE)]
        temp.t <- t(temp)
        colnames(temp.t) <- c("counts")
        temp.df <- as.data.frame(temp.t)
        #k=1
        #OutVals <- boxplot(temp.t,plot=FALSE)$out
        #OutVals <- unique(OutVals)
        #OutVals <- sort(OutVals,decreasing = T)
        #while(k<=length(OutVals)){
        #  # return integer(0), temp.df <- temp.df
        #  if(length(grep(OutVals[k],temp.df$counts))==0){
        #    temp.df <- temp.df
        #  } else {
        #    temp.df <- as.data.frame(temp.df[-grep(OutVals[k],temp.df$counts),])
        #    colnames(temp.df) <- "counts"
        #    #print(length(temp.df$counts))
        #  }
        #  k=k+1
        #}
        temp.df$CancerType<-rep(type[i],times=length(temp.df$counts))
        temp.df$microbe<-rep(target,times=length(temp.df$counts))
        gene <- rbind(gene,temp.df)
        i=i+1
      }
      
      
      #my_comparisons <- list(c("NC_PKU","CRC"),c("NC_PKU","STAD"),c("NC_PKU","LUAD"),c("NC_PKU","ESCA_PKU"),c("NC_ShH","HCC_ShH"),c("NC_ChQ","HCC_ChQ")) #pico
      #my_comparisons <- list(c("NC_PKU","CRC"),c("NC_PKU","STAD"),c("NC_PKU","LUAD"),c("NC_PKU","ESCA_PKU")) #pico_PKU
      #my_comparisons <- list(c("NC_ShH","HCC_ShH")) #pico_ShH
      #my_comparisons <- list(c("NC","CRC"),c("NC","CHD"),c("NC","HCC"),c("NC","PAAD")) #exorbase
      #my_comparisons <- list(c("N","PDAC"),c("N","CP")) #gut
      #my_comparisons <- list(c("HD","NSCLC")) #GSE89843 platelet all batch
      #my_comparisons <- list(c("Vumc.HD","Vumc.NSCLC")) #GSE89843 in vumc
      #my_comparisons <- list(c("HD","Breast"),c("HD","Liver"),c("HD","CRC"),c("HD","Lung"),c("HD","Panc"),c("HD","NSCLC"),c("HD","GBM"))  #GSE68086
      #my_comparisons <- list(c("N","PDAC")) #GSE133684
      #my_comparisons <- list(c("HD","PDAC")) #GSE40174
      #my_comparisons <- list(c("N","C")) #GSE131512
      #my_comparisons <- list(c("N","AD")) #GSE136243 
      #my_comparisons <- list(c("N","HCC")) #GSE120663 
      #my_comparisons <- list(c("NC","CRC")) #PBMC
      #my_comparisons <- list(c("N","THCA")) #THCA
      #my_comparisons <- list(c("NC_PKU","I.CRC"),c("NC_PKU","I.STAD"),c("NC_PKU","I.LUAD"),c("NC_PKU","I.ESCA_PKU"),c("NC_ShH","I.HCC_ShH")) # pico_stage
      #my_comparisons <- list(c("M.NC_PKU","F.NC_PKU"),c("M.CRC","F.CRC"),c("M.STAD","F.STAD"),c("M.LUAD","F.LUAD"),c("M.ESCA_PKU","F.ESCA_PKU"),c("M.HCC_ShH","F.HCC_ShH")) # pico_gender
      #my_comparisons <- list(c("normal","benign"),c("normal","malignant"),c("normal","post.surgery"),c("normal","ectopic")) #GSE27562 PBMC array
      #my_comparisons <- list(c("N","T"))
      #my_comparisons <- list(c("NC","CRC"),c("NC","STAD"))
      #my_comparisons <- list(c("N","T"))
      my_comparisons <- list(c("NC","CRC"))
      
      gene$CancerType <- factor(gene$CancerType,levels=type)
      p <- ggplot(gene,aes(x=CancerType,y=counts,fill=CancerType))+
        geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
        geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
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
          axis.title.y = element_text(face="bold",color="black", size=24))
      
      m=1
      if(m>0){
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  method.args = list(alternative = "greater"),
                                  label = "p.signif"
        )+labs(x="",y="log2(TPM+1)",title=paste0(target,"\nwilcox.test.greater"), face="bold",fill="Type")
        ggsave(p,filename = paste0(gsub("/","",target),"_greater.pdf"),path = "./",width=10,height=8)
      } else if(m==0){
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  method.args = list(alternative = "two.sided"),
                                  label = "p.signif"
        )+labs(x="",y="log2(TPM+1)",title=paste0(target,"\nwilcox.test.twoside"), face="bold",fill="Type")
        ggsave(p,filename = paste0(gsub("/","",target),"_twosided.pdf"),path="./",width=10,height=8)
      } else {
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  method.args = list(alternative = "less"),
                                  label = "p.signif"
        )+labs(x="",y="log2(TPM+1)",title=paste0(target,"\nwilcox.test.less"), face="bold",fill="Type")
        ggsave(p,filename = paste0(gsub("/","",target),"_less.pdf"),path="./",width=10,height=8)
      }
      j=j+1
    }
  }
}


 #rawcounts to TPM
{
  counts <- read.csv("/Users/yuhuan/Desktop/PBMC/GSE120663_all.counts.txt",sep = "\t", header = TRUE, row.names = 1)
  
  annotations <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"), 
                       filters = "hgnc_symbol", 
                       values=rownames(counts), mart= mart,useCache = FALSE)
  library(dplyr)
  longest_transcript <- arrange(annotations, hgnc_symbol, desc(transcript_length))
  
  #merge 会改变顺序，尽量不要用
  #merged <- merge(as.data.frame(rownames(counts)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by.x=1,by.y="hgnc_symbol",all.x=TRUE)
  
  #left_join 不会改变顺序
  merged <- left_join(as.data.frame(rownames(counts)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by = c("rownames(counts)"="hgnc_symbol"))
  
  #对于无注释长度的转录本，记为1000，即不对长度标准化
  merged$transcript_length[is.na(merged$transcript_length)] <- 1000
  
  transcript_lengths <- as.vector(merged$transcript_length)
  
  #find gene length normalized values 
  rpk <- apply( counts, 2, function(x) x/(transcript_lengths/1000))
  #normalize by the sample size using rpk values
  tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
  
  new_rownames <- unite(merged, "new_rownames", `rownames(counts)`, ensembl_gene_id, transcript_length, sep = "|")
  
  rownames(tpm) <- new_rownames$new_rownames
  
  write.csv(tpm,"/Users/yuhuan/Desktop/PBMC/GSE120663_TPM.txt")
}


#dataset comparison heatmap
{
#firstly prepare each Riboosme pathway genes' TPM or log2(TPM+1) matrix based on TPM
{
    #setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/overview/TPM/")
    setwd("/Users/yuhuan/Desktop/Alt/")
    Ribosome <- subset(PATH_ID_NAME,DESCRPTION=="Ribosome")

    target <- unique(Ribosome$ensembl_gene_id)
    #pico
    #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/TPM_all_pico.txt",sep="\t",header = T, row.names = 1)
    #type <- c("NC_ShH","HCC_ShH")
    #type <- c("NC_ChQ","HCC_ChQ")
    #type <- c("NC_PKU","CRC","STAD","LUAD","ESCA_PKU")
    #exoRbase
    #type <- c("NC","CRC","CHD","HCC","PAAD")
    #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/exoRbase_TPM.txt",sep="\t",header = T, row.names = 1)
    #GSE133684
    #type <- c("N","CP","PDAC")
    #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/GSE133684_exp_TPM-all.txt",sep="\t",header = T, row.names = 1)
    #GSE89843 vumc
    #type <- c("Vumc.HD","Vumc.NSCLC")
    #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/GSE89843_TPM.txt",sep="\t",header = T, row.names = 1)
    #GSE68086
    #type <- c("HD","Breast","Liver","CRC","Lung","Panc","NSCLC","GBM")
    #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/GSE68086_TPM.txt",sep="\t",header = T, row.names = 1)
    #GSE117623
    #type <- c("CLD","HCC")
    #counts <- read.csv("/Users/yuhuan/Desktop/CEC/GSE117623_TPM.csv",header = T, row.names = 1)
    #GSE144561
    #type <- c("HD","PDAC")
    #counts <- read.csv("/Users/yuhuan/Desktop/CTC/GSE144561_TPM.csv",header = T, row.names = 1)
    #GSE106804
    #type <- c("HD","GBM")
    #counts <- read.csv("/Users/yuhuan/Desktop/EV/GSE106804_TPM.csv",header = T, row.names = 1)
    #GSE120663
    #type <- c("N","HCC")
    #counts <- read.csv("/Users/yuhuan/Desktop/PBMC/GSE120663_TPM.csv",header = T, row.names = 1)
    #pico_PBMC(CRC)
    #type <- c("N","T","PBMC","plasma","HD")
    #counts <- read.csv("/Users/yuhuan/Desktop/PBMC/pico_PBMC_TPM.csv",header = T, row.names = 1)
    #GSE133684_promoter
    #type <- c("SRR")
    #counts <- read.csv("/Users/yuhuan/Desktop/Alt/TPM-by-promoter_GSE133684.txt",sep = "\t",header = T, row.names = 1)
    
    
    #counts <- log2(as.matrix(counts)+1)
    
    df.type=as.data.frame(array(dim=c(1,1))) #设置数据框行列
    i=1
    while(i<=length(type)){
      df <- counts[,grep(type[i],colnames(counts),fixed=TRUE)]
      df.type <- cbind(df.type,df)
      i=i+1
    }
    df.type <- df.type[,-1]
    
    # progress
    pb <- progress_bar$new(
      format = "  Processing [:bar] :percent eta: :eta",
      total = length(target), clear = FALSE, width= 60) 
    df.target={}
    j=1
    while(j<=length(target)){
      df <- df.type[grep(target[j],rownames(df.type),fixed=TRUE),]
      if(nrow(df)==0){
        df <- data.frame(matrix(nrow=1,ncol=ncol(df.type)))
        rownames(df) <- target[j]
        colnames(df) <- colnames(df.type)
        df.target <- rbind(df.target,df)
      }else{
        df.target <- rbind(df.target,df)
      }
      j=j+1
      pb$tick()
      Sys.sleep(1 / 100)
    }

write.csv(df.target,"TPM-by-promoter_GSE133684_ribosome.csv")

}

#manually summarize and draw heatmap
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/overview/log2/")

target <- unique(Ribosome$ensembl_gene_id)
  
ribosome_datasets <- list.files('./',pattern="*.csv")

readin_datasets <- lapply(ribosome_datasets, function(x) read.csv(x, header = TRUE, row.names = 1)) 

overview <- do.call("cbind", readin_datasets)

genes <- merge(as.data.frame(target),Ribosome[,-which(colnames(Ribosome)=="ENTREZID")],by.x="target",by.y="ensembl_gene_id",all.x=TRUE)

genes <- genes[!duplicated(genes),]

genes <- unite(genes,"rownames",target,hgnc_symbol,remove=FALSE,sep = ":")

rownames(overview) <- genes$rownames

write.csv(overview,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/overview/overview_log2.csv")
}
  
{
overview <- read.csv("./heatmap/overview_log2.csv",header = TRUE, row.names = 1)
col_annotation <- read.csv("./heatmap/col_annotation_log2.csv",header = T, row.names = 1 )
ann_colors = list(
  Groups = c(HD = "green", pre_Ca = "yellow", CP = "#5CACEE",PAAD = "#4E78A0", PDAC = "#4E78A0", CRC = "#0D4F8B", CHD = "#104E8B", HCC = "#003F87", STAD = "#2B4F81", LUAD = "#191970", ESCA = "#27408B", NSCLC = "#191970", GBM = "#3D59AB", BRCA = "#26466D"),
  #Datatsets = c(exoRBase = "#ECC3BF", GSE133684 = "#C75D4D", GSE68086 = "#FF3030", GSE89843 = "#FC1501", pico_PKU = "#F4A460",pico_ShH = "#FEE8D6",pico_ChQ = "#FCD59C"),
  Specimen = c(EV = "#FFCCCC", plasma = "#FFA824", platelet = "#CC1100",adjacent_tissue="blue",tumor="red",PBMC="white",CEC="pink",CTC="purple")
)


bk = unique(c(seq(-0.5,0.5, length=100)))
pheatmap(overview,
         scale = "row",
         annotation_col = col_annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         breaks=bk,
         na_col = "white",
         colorRampPalette(c("blue","white","red"))(100),
         show_colnames=FALSE, 
         fontsize_row = 8,gaps_col = c(30,168,249,261,882,1487))

#ggsave(p,"test.pdf",width=50,height=30)
}
}


##coefficent of variation
{
  library(ineq)
  TPM <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/TPM_all_pico.txt",sep = "\t",header = T, row.names = 1)
  
  #cv
  sd <- rowSdDiffs(TPM)
  mean <- rowMeans(TPM)
  cv <- sd/mean
  
  #gini by ineq
  gini <- apply(TPM, MARGIN = 1, function(x) Gini(x))
  
  #cv by ineq
  test <- apply(TPM, MARGIN = 1, function(x) var.coeff(x))
  
  #plot cv2mean
  #cv2mean <- cbind(as.data.frame(mean),as.data.frame(cv))
  #ggplot(cv2mean,aes(x=mean,y=cv))+geom_point()+theme_bw()+xlab("Average TPM")+ylab("Coefficient of Variation")+ylim(0,5)+xlim(0,50000)
  
  #plot gini2mean
  gini2mean <- cbind(as.data.frame(mean),as.data.frame(gini))
  
  gini2mean$symbol <- as.character(lapply(strsplit(rownames(gini2mean),"\\|"), function(x) x[3]))
  label <- as.data.frame(Ribosome_symbol,Ribosome_symbol)
  label$color <- "Ribosome"
  plot <- left_join(gini2mean,label,by=c("symbol"="Ribosome_symbol"))
  plot$color[is.na(plot$color)] <- "Other"
  top50 <- head(plot[order(plot$mean,decreasing = TRUE),c(1,3)],50)
  colnames(top50) <- c("mean","label")
  plot <- left_join(plot,top50,by=c("mean"="mean"))
  
  plot[which(plot$symbol=="GAPDH"),c("label")] <- "GAPDH"
  
  plot[which(plot$symbol=="GAPDH"),c("color")] <- "Inter Control"
  plot[which(plot$symbol=="ACTB"),c("color")] <- "Inter Control"
  
  ggplot(plot,aes(x=mean,y=gini,color=color))+
    scale_color_manual(values = c("blue","grey","red"))+
    geom_point()+
    theme_bw()+
    xlab("Average TPM")+
    ylab("gini")+
    geom_text_repel(aes(x=mean,y=gini,label=label))
}

#Alt promoter heatmap
{
  overview <- read.csv("./TPM-by-promoter_GSE133684_ribosome.csv",header = TRUE, row.names = 1)
  col_annotation <- read.csv("./col_annotation_promoter_GSE133684.csv",header = T, row.names = 1 )
  #ann_colors = list(
  #  Groups = c(HD = "green", pre_Ca = "yellow", CP = "#5CACEE",PAAD = "#4E78A0", PDAC = "#4E78A0", CRC = "#0D4F8B", CHD = "#104E8B", HCC = "#003F87", STAD = "#2B4F81", LUAD = "#191970", ESCA = "#27408B", NSCLC = "#191970", GBM = "#3D59AB", BRCA = "#26466D"),
  #  #Datatsets = c(exoRBase = "#ECC3BF", GSE133684 = "#C75D4D", GSE68086 = "#FF3030", GSE89843 = "#FC1501", pico_PKU = "#F4A460",pico_ShH = "#FEE8D6",pico_ChQ = "#FCD59C"),
  #  Specimen = c(EV = "#FFCCCC", plasma = "#FFA824", platelet = "#CC1100",adjacent_tissue="blue",tumor="red",PBMC="white",CEC="pink",CTC="purple")
  #)
  bk = unique(c(seq(-0.5,0.5, length=100)))
  pheatmap(overview,
           scale = "row",
           annotation_col = col_annotation,
           #annotation_colors = ann_colors,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           breaks=bk,
           na_col = "white",
           colorRampPalette(c("blue","white","red"))(100),
           show_colnames=FALSE, 
           fontsize_row = 8
           )
}

#plot MT-RNA
{
  counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/TPM_all_pico.txt",sep="\t",header = T, row.names = 1)
  counts <- log2(as.matrix(counts)+1)
  MT={}
  MT$ensembl_gene_id <- c("ENSG00000198888","ENSG00000198763","ENSG00000198840","ENSG00000212907","ENSG00000198886",
                           "ENSG00000198786","ENSG00000198695","ENSG00000198727","ENSG00000198804","ENSG00000198712",
                           "ENSG00000198938","ENSG00000198899","ENSG00000228253","ENSG00000211459","ENSG00000210082",
                           "ENSG00000210127","ENSG00000210174","ENSG00000210135","ENSG00000210154","ENSG00000210140",
                           "ENSG00000210194","ENSG00000210107","ENSG00000210164","ENSG00000210176","ENSG00000210100",
                           "ENSG00000209082","ENSG00000210191","ENSG00000210156","ENSG00000210112","ENSG00000210049",
                           "ENSG00000210196","ENSG00000210151","ENSG00000210184","ENSG00000210195","ENSG00000210117",
                           "ENSG00000210144","ENSG00000210077")
  MT$gene_symbol <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-ND6","MT-CYB","MT-CO1","MT-CO2","MT-CO3",
                      "MT-ATP6","MT-ATP8","MT-RNR1","MT-RNR2","MT-TA","MT-TR","MT-TN","MT-TD","MT-TC","MT-TE","MT-TQ","MT-TG","MT-TH",
                      "MT-TI","MT-TL1","MT-TL2","MT-TK","MT-TM","MT-TF","MT-TP","MT-TS1","MT-TS2","MT-TT","MT-TW","MT-TY","MT-TV")
  
  MT <- as.data.frame(MT)
  
  MT$DESCRPTION <- "Mitochondrial RNA"
  
  # set progress
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = length(unique(MT$DESCRPTION)), clear = FALSE, width= 60) 
  #get all genes in each pathway in KEGG
  n=1
  pathway_count={}
  while(n <= length(unique(MT$DESCRPTION))){
    pathway_name <- as.character(unique(MT$DESCRPTION)[n])
    pathway <- subset(MT,DESCRPTION==pathway_name)
    
    candidate <- pathway$ensembl_gene_id #for colnames with ensemble id
    #candidate <- pathway$hgnc_symbol #for colnames with gene symbol
    
    candidate <- unique(candidate)
    
    #summary one pathway total count in each sample and make 1 row dataframe
    j=1
    pathway_gene_count={}
    while(j<=length(candidate)){
      target <- candidate[j]
      if(length(grep(target,rownames(counts)))==0) {
        #print(paste0("No ",target," in this dataset."))
        j=j+1
      } else {
        #temp <- counts[which(rownames(counts)==target),]
        temp <- counts[grep(target,rownames(counts),fixed=TRUE),]
        pathway_gene_count <- rbind(pathway_gene_count,temp)
        j=j+1
      }
    }
    if(is.null(pathway_gene_count)){
      n=n+1
      next;
    } else {
      one_pathway <- as.data.frame(t(colMeans(pathway_gene_count))) #colSums
      rownames(one_pathway) <- pathway_name
      pathway_count <- rbind(pathway_count,one_pathway)
      n=n+1
    }
    pb$tick()
    Sys.sleep(1 / 100)
    
  }
  write.csv(pathway_count,"pico_mtRNA_TPM_pathway_level_log2.csv")
}

#plot top genes in RNAseq
{
#preparation
{
  library(ggplot2)
  library(tidyverse)
  library(extrafont)
  fonts()
}

{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/MT-RNA")
  counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/TPM_all_pico.txt",sep="\t",header = T, row.names = 1)
  #phylum <- read.csv("genus2phylum.csv",header= TRUE,row.names = 1)
  type <- c("NC","CRC","STAD","HCC","LUAD","ESCA")
  top=50
  #color=c("Actinobacteria"="red","Tenericutes"="yellow","Bacteroidetes"="blue",
  #        "Cyanobacteria"="green","Proteobacteria"="purple","Firmicutes"="brown",
  #        "Virus"="grey","Unknown"="black","Euryarchaeota"="white")
}


j=1
while(j<=length(type)){
  counts_type <- counts[,grep(type[j],colnames(counts),fixed=TRUE)]
  counts_type$sum <- rowSums(counts_type)
  counts_type.order <- counts_type[order(counts_type$sum,decreasing=T),]
  
  candidate <- rownames(counts_type.order)
  i=1
  microbe={}
  while(i<=top){
    temp <- counts_type.order[grep(candidate[i],rownames(counts_type.order),fixed=TRUE)[1],!names(counts_type.order) %in% c("sum")]
    temp.t <- t(temp)
    colnames(temp.t) <- c("counts")
    temp.df <- as.data.frame(temp.t)
    #k=1
    #OutVals <- boxplot(temp.t,plot=FALSE)$out
    #OutVals <- unique(OutVals)
    #OutVals <- sort(OutVals,decreasing = T)
    #while(k<=length(OutVals)){
    #  # return integer(0), temp.df <- temp.df
    #  if(length(grep(OutVals[k],temp.df$counts))==0){
    #    temp.df <- temp.df
    #  } else {
    #    temp.df <- as.data.frame(temp.df[-grep(OutVals[k],temp.df$counts),])
    #    colnames(temp.df) <- "counts"
    #    #print(length(temp.df$counts))
    #  }
    #  k=k+1
    #}
    temp.df$microbe<-rep(candidate[i],times=length(temp.df$counts))
    #temp.df$phylum<-rep(phylum[grep(candidate[i],rownames(phylum),fixed=TRUE),],times=length(temp.df$counts))
    microbe <- rbind(microbe,temp.df)
    i=i+1
  }
  
  level <- as.character(lapply(strsplit(candidate[1:top],"|",fixed = TRUE), function(x) x[3]))
  microbe$symbol <- as.character(lapply(strsplit(as.character(microbe$microbe),"|",fixed=TRUE), function(x) x[3]))
  microbe$symbol <- factor(microbe$symbol,levels = level)
  p <- ggplot(microbe,aes(x=symbol,y=counts))+
    geom_boxplot(alpha = 1, size=0, position = position_dodge(1.1),outlier.size=2,outlier.alpha = 0.6,outlier.color = "grey")+
    #geom_point(size = 2, color="grey", alpha = 0.6)+
    #scale_fill_manual(values=color,name="Phylum")+
    #scale_fill_brewer(palette="Blues") +
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = margin(2,2,2,50,unit="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1,vjust = 1),
      #axis.ticks.x = element_line(hjust=0),
      axis.ticks.length.x = unit(-0.2,"cm"),
      axis.text.y = element_text(color="black", size=24),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    labs(x="",y="TPM",title=type[j])+
    stat_summary(fun.y=mean, geom="line",linetype="dotted",
                 size=1, color = "black", alpha =0.6, aes(group=1))+
    stat_summary(fun.y=mean, geom="point",alpha=0.8, color = "black",
                 size = 0.5, shape= 21, stroke = 1)
  
  p
  ggsave(p,filename=paste0(type[j],"_Top",top,".pdf"),path="./")
  j=j+1
}
}

# percent of MT-RNA(sumed TPM/ 1M) in RNAseq
{
  counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/TPM_all_pico.txt",sep="\t",header = T, row.names = 1)
  #counts <- log2(as.matrix(counts)+1)
  MT={}
  MT$ensembl_gene_id <- c("ENSG00000198888","ENSG00000198763","ENSG00000198840","ENSG00000212907","ENSG00000198886",
                          "ENSG00000198786","ENSG00000198695","ENSG00000198727","ENSG00000198804","ENSG00000198712",
                          "ENSG00000198938","ENSG00000198899","ENSG00000228253","ENSG00000211459","ENSG00000210082",
                          "ENSG00000210127","ENSG00000210174","ENSG00000210135","ENSG00000210154","ENSG00000210140",
                          "ENSG00000210194","ENSG00000210107","ENSG00000210164","ENSG00000210176","ENSG00000210100",
                          "ENSG00000209082","ENSG00000210191","ENSG00000210156","ENSG00000210112","ENSG00000210049",
                          "ENSG00000210196","ENSG00000210151","ENSG00000210184","ENSG00000210195","ENSG00000210117",
                          "ENSG00000210144","ENSG00000210077")
  MT$gene_symbol <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-ND6","MT-CYB","MT-CO1","MT-CO2","MT-CO3",
                      "MT-ATP6","MT-ATP8","MT-RNR1","MT-RNR2","MT-TA","MT-TR","MT-TN","MT-TD","MT-TC","MT-TE","MT-TQ","MT-TG","MT-TH",
                      "MT-TI","MT-TL1","MT-TL2","MT-TK","MT-TM","MT-TF","MT-TP","MT-TS1","MT-TS2","MT-TT","MT-TW","MT-TY","MT-TV")
  
  MT <- as.data.frame(MT)
  
  MT$DESCRPTION <- "Mitochondrial RNA"
  
  # set progress
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = length(unique(MT$DESCRPTION)), clear = FALSE, width= 60) 
  #get all genes in each pathway in KEGG
  n=1
  pathway_count={}
  while(n <= length(unique(MT$DESCRPTION))){
    pathway_name <- as.character(unique(MT$DESCRPTION)[n])
    pathway <- subset(MT,DESCRPTION==pathway_name)
    
    candidate <- pathway$ensembl_gene_id #for colnames with ensemble id
    #candidate <- pathway$hgnc_symbol #for colnames with gene symbol
    
    candidate <- unique(candidate)
    
    #summary one pathway total count in each sample and make 1 row dataframe
    j=1
    pathway_gene_count={}
    while(j<=length(candidate)){
      target <- candidate[j]
      if(length(grep(target,rownames(counts)))==0) {
        #print(paste0("No ",target," in this dataset."))
        j=j+1
      } else {
        #temp <- counts[which(rownames(counts)==target),]
        temp <- counts[grep(target,rownames(counts),fixed=TRUE),]
        pathway_gene_count <- rbind(pathway_gene_count,temp)
        j=j+1
      }
    }
    if(is.null(pathway_gene_count)){
      n=n+1
      next;
    } else {
      one_pathway <- as.data.frame(t(colSums(pathway_gene_count)/1000000)) #colSums
      rownames(one_pathway) <- pathway_name
      pathway_count <- rbind(pathway_count,one_pathway)
      n=n+1
    }
    pb$tick()
    Sys.sleep(1 / 100)
    
  }
  write.csv(pathway_count,"pico_mtRNA_TPM_pathway_level_percent_in_RNAseq.csv")
}

#clustering by 78 RP-mRNA in NAR
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/clustering/")
  #pico
  #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/Raw_counts_passed_QC.txt",sep="\t",header = T, row.names = 1)
  #GSE68086
  #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/GSE68086_TEP_data_matrix_forDE.txt",sep="\t",header = T, row.names = 1)
  #PBMC
  #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/PBMC.txt",sep="\t",header = T, row.names = 1)
  #exoRBase
  #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/exoRbase_rawcounts.txt",sep="\t",header = T, row.names = 1)
  #GSE133684
  counts <- read.csv("/Users/yuhuan/Desktop/EV/GSE133684.txt",sep="\t",header = T, row.names = 1)
  #pico_tissue
  #counts <- read.csv("/Users/yuhuan/Desktop/tissue/pico_tissue_rawcount.txt",sep="\t",header = T, row.names = 1)
  
  #not success
  #PBMC GSE120663
  #counts <- read.csv("/Users/yuhuan/Desktop/PBMC/GSE120663_all.counts.txt",sep="\t",header = T, row.names = 1)
  #CTC GSE144561
  #counts <- read.csv("/Users/yuhuan/Desktop/CTC/GSE144561_rawCountsAllsamples.txt",sep="\t",header = T, row.names = 1)
  #CEC GSE117623
  #counts <- read.csv("/Users/yuhuan/Desktop/CEC/GSE117623_RawCounts_Geo_CLDHCCWBC.csv",header = T, row.names = 1)
  RP_RNA <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/RP-RNA.csv",header = T)
  
  #get 78 RP RNAs
  j=1
  pathway_gene_count={}
  while(j<=nrow(RP_RNA)){
    target <- RP_RNA[j,2]
    gene_symbol <- RP_RNA[j,1]
    if(length(grep(target,rownames(counts)))==0) {
      #print(paste0("No ",target," in this dataset."))
      temp <- as.data.frame(array(,dim=c(1,ncol(counts))))
      temp[1,] <- 0
      rownames(temp) <- gene_symbol
      j=j+1
    } else {
      #temp <- counts[which(rownames(counts)==target),]
      #temp <- counts[grep(target,rownames(counts),fixed=TRUE),]  #for ensg
      temp <- counts[grep(target,rownames(counts),fixed=TRUE),] #for gene symbol
      rownames(temp) <- gene_symbol
      pathway_gene_count <- rbind(pathway_gene_count,temp)
      j=j+1
    }
  }
  write.csv(pathway_gene_count,"GSE133684_RP-RNA.csv")
  
  #normalize
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  annotations <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"), 
                       filters = "hgnc_symbol", 
                       values=rownames(pathway_gene_count), mart= mart,useCache = FALSE)
  library(dplyr)
  longest_transcript <- arrange(annotations, hgnc_symbol, desc(transcript_length))
  
  #merge 会改变顺序，尽量不要用
  #merged <- merge(as.data.frame(rownames(counts)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by.x=1,by.y="hgnc_symbol",all.x=TRUE)
  
  #read local RP-RNA count
  pathway_gene_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/clustering/different_component.csv",header = TRUE, row.names = 1)
  #left_join 不会改变顺序
  merged <- left_join(as.data.frame(rownames(pathway_gene_count)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by = c("rownames(pathway_gene_count)"="hgnc_symbol"))
  
  #对于无注释长度的转录本，记为1000，即不对长度标准化
  merged$transcript_length[is.na(merged$transcript_length)] <- 1000
  
  transcript_lengths <- as.vector(merged$transcript_length)
  
  #find gene length normalized values 
  rpk <- apply(pathway_gene_count, 2, function(x) x/(transcript_lengths/1000))
  #normalize by the sample size using rpk values
  normalized <- apply(rpk, 2, function(x) x / sum(as.numeric(x)))
  
  new_rownames <- unite(merged, "new_rownames", `rownames(pathway_gene_count)`, ensembl_gene_id, transcript_length, sep = "|")
  
  rownames(normalized) <- new_rownames$new_rownames
  
  #t-sne
  library(Rtsne)
  
  normalized.t <- as.data.frame(t(normalized))
  ##nrow(normlized.t) larger than 3*perplexity
  tsne_out <- Rtsne(normalized.t,pca=FALSE,dims=2,perplexity = 30,theta=0.0) 
  
  tsne_res <- as.data.frame(tsne_out$Y)
  colnames(tsne_res) <- c("tSNE1","tSNE2")
  rownames(tsne_res) <- rownames(normalized.t)
  tsne_res$group <- as.character(lapply(strsplit(rownames(normalized.t),".",fixed = TRUE),function(x) x[1]))
  #tsne_res$group <- as.character(lapply(strsplit(tsne_res$group,"_",fixed = TRUE),function(x) x[1]))
  head(tsne_res)
  
  ## 拟杆菌Bacteroides、普氏菌Prevotella和瘤胃球菌Ruminococcus
  
  ggplot(tsne_res,aes(tSNE1,tSNE2,color=tsne_res$group)) + 
    geom_point() + theme_bw() + 
    #scale_color_gradient(low = "blue",high = "red",) +
    #geom_hline(yintercept = 0,lty=2,col="red") + 
    #geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
    #stat_ellipse(aes(tSNE1,tSNE2, color = tsne_res$group), geom = "polygon", alpha = 0.3, levels = 0.99)+
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(title = "RP-RNA")
  
}

#clustering by 78 RP-mRNA in NAR for TCGA
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/clustering/TCGA/")
  counts <- read.csv("THCA.txt",sep="\t",header = T, row.names = 1)
  RP_RNA <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/RP-RNA.csv",header = T)
  
  #get 78 RP RNAs
  j=1
  pathway_gene_count={}
  while(j<=nrow(RP_RNA)){
    target <- RP_RNA[j,2]
    gene_symbol <- RP_RNA[j,1]
    if(length(grep(target,rownames(counts)))==0) {
      #print(paste0("No ",target," in this dataset."))
      temp <- as.data.frame(array(,dim=c(1,ncol(counts))))
      temp[1,] <- 0
      rownames(temp) <- gene_symbol
      j=j+1
    } else {
      #temp <- counts[which(rownames(counts)==target),]
      temp <- counts[grep(target,rownames(counts),fixed=TRUE),]
      rownames(temp) <- gene_symbol
      pathway_gene_count <- rbind(pathway_gene_count,temp)
      j=j+1
    }
  }
  write.csv(pathway_gene_count,"THCA_RP-RNA.csv")
  
  pathway_gene_count <- read.csv("TCGA_RP.csv",header=TRUE,row.names = 1)
  #normalize
  #mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  #annotations <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"), 
  #                     filters = "hgnc_symbol", 
  #                     values=rownames(pathway_gene_count), mart= mart,useCache = FALSE)
  library(dplyr)
  longest_transcript <- arrange(annotations, hgnc_symbol, desc(transcript_length))
  
  #merge 会改变顺序，尽量不要用
  #merged <- merge(as.data.frame(rownames(counts)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by.x=1,by.y="hgnc_symbol",all.x=TRUE)
  
  #left_join 不会改变顺序
  merged <- left_join(as.data.frame(rownames(pathway_gene_count)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by = c("rownames(pathway_gene_count)"="hgnc_symbol"))
  
  #对于无注释长度的转录本，记为1000，即不对长度标准化
  merged$transcript_length[is.na(merged$transcript_length)] <- 1000
  
  transcript_lengths <- as.vector(merged$transcript_length)
  
  #find gene length normalized values 
  rpk <- apply(pathway_gene_count, 2, function(x) x/(transcript_lengths/1000))
  #normalize by the sample size using rpk values
  normalized <- apply(rpk, 2, function(x) x / sum(as.numeric(x)))
  
  new_rownames <- unite(merged, "new_rownames", `rownames(pathway_gene_count)`, ensembl_gene_id, transcript_length, sep = "|")
  
  rownames(normalized) <- new_rownames$new_rownames
  
  #t-sne
  library(Rtsne)
  
  normalized.t <- as.data.frame(t(normalized))
  ##nrow(normlized.t) larger than 3*perplexity
  tsne_out <- Rtsne(normalized.t,pca=FALSE,dims=2,perplexity = 300,theta=0.0) 
  
  tsne_res <- as.data.frame(tsne_out$Y)
  colnames(tsne_res) <- c("tSNE1","tSNE2")
  rownames(tsne_res) <- rownames(normalized.t)
  tsne_res$group <- as.character(lapply(strsplit(rownames(normalized.t),"_",fixed = TRUE),function(x) x[1]))
  #tsne_res$group <- as.character(lapply(strsplit(tsne_res$group,"_",fixed = TRUE),function(x) x[1]))
  head(tsne_res)

  ggplot(tsne_res,aes(tSNE1,tSNE2,color=tsne_res$group)) + 
    geom_point() + theme_bw() + 
    #scale_color_gradient(low = "blue",high = "red",) +
    #geom_hline(yintercept = 0,lty=2,col="red") + 
    #geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
    #stat_ellipse(aes(tSNE1,tSNE2, color = tsne_res$group), geom = "polygon", alpha = 0.3, levels = 0.99)+
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(title = "RP-RNA")
  
}

#clustering by 37 MT RNA
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/clustering/")
  counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/Raw_counts_passed_QC.txt",sep="\t",header = T, row.names = 1)
  MT={}
  MT$ensembl_gene_id <- c("ENSG00000198888","ENSG00000198763","ENSG00000198840","ENSG00000212907","ENSG00000198886",
                          "ENSG00000198786","ENSG00000198695","ENSG00000198727","ENSG00000198804","ENSG00000198712",
                          "ENSG00000198938","ENSG00000198899","ENSG00000228253","ENSG00000211459","ENSG00000210082",
                          "ENSG00000210127","ENSG00000210174","ENSG00000210135","ENSG00000210154","ENSG00000210140",
                          "ENSG00000210194","ENSG00000210107","ENSG00000210164","ENSG00000210176","ENSG00000210100",
                          "ENSG00000209082","ENSG00000210191","ENSG00000210156","ENSG00000210112","ENSG00000210049",
                          "ENSG00000210196","ENSG00000210151","ENSG00000210184","ENSG00000210195","ENSG00000210117",
                          "ENSG00000210144","ENSG00000210077")
  MT$gene_symbol <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-ND6","MT-CYB","MT-CO1","MT-CO2","MT-CO3",
                      "MT-ATP6","MT-ATP8","MT-RNR1","MT-RNR2","MT-TA","MT-TR","MT-TN","MT-TD","MT-TC","MT-TE","MT-TQ","MT-TG","MT-TH",
                      "MT-TI","MT-TL1","MT-TL2","MT-TK","MT-TM","MT-TF","MT-TP","MT-TS1","MT-TS2","MT-TT","MT-TW","MT-TY","MT-TV")
  
  MT <- as.data.frame(MT)
  
  MT$DESCRPTION <- "Mitochondrial RNA"
  
  #get MT RNAs
  j=1
  pathway_gene_count={}
  while(j<=nrow(MT)){
    target <- MT[j,1]
    gene_symbol <- MT[j,2]
    if(length(grep(target,rownames(counts)))==0) {
      #print(paste0("No ",target," in this dataset."))
      temp <- as.data.frame(array(,dim=c(1,ncol(counts))))
      temp[1,] <- 0
      rownames(temp) <- gene_symbol
      j=j+1
    } else {
      #temp <- counts[which(rownames(counts)==target),]
      temp <- counts[grep(target,rownames(counts),fixed=TRUE),]
      rownames(temp) <- gene_symbol
      pathway_gene_count <- rbind(pathway_gene_count,temp)
      j=j+1
    }
  }
  write.csv(pathway_gene_count,"pico_MT-RNA.csv")
  
  #normalize
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  annotations <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"), 
                       filters = "hgnc_symbol", 
                       values=rownames(pathway_gene_count), mart= mart,useCache = FALSE)
  library(dplyr)
  longest_transcript <- arrange(annotations, hgnc_symbol, desc(transcript_length))
  
  #merge 会改变顺序，尽量不要用
  #merged <- merge(as.data.frame(rownames(counts)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by.x=1,by.y="hgnc_symbol",all.x=TRUE)
  
  #left_join 不会改变顺序
  merged <- left_join(as.data.frame(rownames(pathway_gene_count)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by = c("rownames(pathway_gene_count)"="hgnc_symbol"))
  
  #对于无注释长度的转录本，记为1000，即不对长度标准化
  merged$transcript_length[is.na(merged$transcript_length)] <- 1000
  
  transcript_lengths <- as.vector(merged$transcript_length)
  
  #find gene length normalized values 
  rpk <- apply(pathway_gene_count, 2, function(x) x/(transcript_lengths/1000))
  #normalize by the sample size using rpk values
  normalized <- apply(rpk, 2, function(x) x / sum(as.numeric(x)))
  
  new_rownames <- unite(merged, "new_rownames", `rownames(pathway_gene_count)`, ensembl_gene_id, transcript_length, sep = "|")
  
  rownames(normalized) <- new_rownames$new_rownames
  
  #t-sne
  library(Rtsne)
  
  normalized.t <- as.data.frame(t(normalized))
  ##nrow(normlized.t) larger than 3*perplexity
  tsne_out <- Rtsne(normalized.t,pca=FALSE,dims=2,perplexity = 120,theta=0.0) 
  
  tsne_res <- as.data.frame(tsne_out$Y)
  colnames(tsne_res) <- c("tSNE1","tSNE2")
  rownames(tsne_res) <- rownames(normalized.t)
  tsne_res$group <- as.character(lapply(strsplit(rownames(normalized.t),".",fixed = TRUE),function(x) x[1]))
  tsne_res$group <- as.character(lapply(strsplit(tsne_res$group,"_",fixed = TRUE),function(x) x[1]))
  head(tsne_res)
  
  ## 拟杆菌Bacteroides、普氏菌Prevotella和瘤胃球菌Ruminococcus
  
  ggplot(tsne_res,aes(tSNE1,tSNE2,color=tsne_res$group)) + 
    geom_point() + theme_bw() + 
    #scale_color_gradient(low = "blue",high = "red",) +
    #geom_hline(yintercept = 0,lty=2,col="red") + 
    #geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
    #stat_ellipse(aes(tSNE1,tSNE2, color = tsne_res$group), geom = "polygon", alpha = 0.3, levels = 0.99)+
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(title = "MT-RNA")
  
}

#clustering by 37 MT RNA for TCGA
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/clustering/TCGA/")
  counts <- read.csv("STAD.txt",sep="\t",header = T, row.names = 1)
  MT={}
  MT$ensembl_gene_id <- c("ENSG00000198888","ENSG00000198763","ENSG00000198840","ENSG00000212907","ENSG00000198886",
                          "ENSG00000198786","ENSG00000198695","ENSG00000198727","ENSG00000198804","ENSG00000198712",
                          "ENSG00000198938","ENSG00000198899","ENSG00000228253","ENSG00000211459","ENSG00000210082",
                          "ENSG00000210127","ENSG00000210174","ENSG00000210135","ENSG00000210154","ENSG00000210140",
                          "ENSG00000210194","ENSG00000210107","ENSG00000210164","ENSG00000210176","ENSG00000210100",
                          "ENSG00000209082","ENSG00000210191","ENSG00000210156","ENSG00000210112","ENSG00000210049",
                          "ENSG00000210196","ENSG00000210151","ENSG00000210184","ENSG00000210195","ENSG00000210117",
                          "ENSG00000210144","ENSG00000210077")
  MT$gene_symbol <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-ND6","MT-CYB","MT-CO1","MT-CO2","MT-CO3",
                      "MT-ATP6","MT-ATP8","MT-RNR1","MT-RNR2","MT-TA","MT-TR","MT-TN","MT-TD","MT-TC","MT-TE","MT-TQ","MT-TG","MT-TH",
                      "MT-TI","MT-TL1","MT-TL2","MT-TK","MT-TM","MT-TF","MT-TP","MT-TS1","MT-TS2","MT-TT","MT-TW","MT-TY","MT-TV")
  
  MT <- as.data.frame(MT)
  
  MT$DESCRPTION <- "Mitochondrial RNA"
  
  #get MT RNAs
  j=1
  pathway_gene_count={}
  while(j<=nrow(MT)){
    target <- MT[j,1]
    gene_symbol <- MT[j,2]
    if(length(grep(target,rownames(counts)))==0) {
      #print(paste0("No ",target," in this dataset."))
      temp <- as.data.frame(array(,dim=c(1,ncol(counts))))
      temp[1,] <- 0
      rownames(temp) <- gene_symbol
      j=j+1
    } else {
      #temp <- counts[which(rownames(counts)==target),]
      temp <- counts[grep(target,rownames(counts),fixed=TRUE),]
      rownames(temp) <- gene_symbol
      pathway_gene_count <- rbind(pathway_gene_count,temp)
      j=j+1
    }
  }
  write.csv(pathway_gene_count,"STAD_MT-RNA.csv")
  
  pathway_gene_count <- read.csv("TCGA_MT.csv",header=TRUE,row.names = 1)
  #normalize
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  annotations <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"), 
                       filters = "hgnc_symbol", 
                       values=rownames(pathway_gene_count), mart= mart,useCache = FALSE)
  library(dplyr)
  longest_transcript <- arrange(annotations, hgnc_symbol, desc(transcript_length))
  
  #merge 会改变顺序，尽量不要用
  #merged <- merge(as.data.frame(rownames(counts)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by.x=1,by.y="hgnc_symbol",all.x=TRUE)
  
  #left_join 不会改变顺序
  merged <- left_join(as.data.frame(rownames(pathway_gene_count)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by = c("rownames(pathway_gene_count)"="hgnc_symbol"))
  
  #对于无注释长度的转录本，记为1000，即不对长度标准化
  merged$transcript_length[is.na(merged$transcript_length)] <- 1000
  
  transcript_lengths <- as.vector(merged$transcript_length)
  
  #find gene length normalized values 
  rpk <- apply(pathway_gene_count, 2, function(x) x/(transcript_lengths/1000))
  #normalize by the sample size using rpk values
  normalized <- apply(rpk, 2, function(x) x / sum(as.numeric(x)))
  
  new_rownames <- unite(merged, "new_rownames", `rownames(pathway_gene_count)`, ensembl_gene_id, transcript_length, sep = "|")
  
  rownames(normalized) <- new_rownames$new_rownames
  
  #t-sne
  library(Rtsne)
  
  normalized.t <- as.data.frame(t(normalized))
  ##nrow(normlized.t) larger than 3*perplexity
  tsne_out <- Rtsne(normalized.t,pca=FALSE,dims=2,perplexity = 30,theta=0.0) 
  
  tsne_res <- as.data.frame(tsne_out$Y)
  colnames(tsne_res) <- c("tSNE1","tSNE2")
  rownames(tsne_res) <- rownames(normalized.t)
  tsne_res$group <- as.character(lapply(strsplit(rownames(normalized.t),"_",fixed = TRUE),function(x) x[1]))
  #tsne_res$group <- as.character(lapply(strsplit(tsne_res$group,"_",fixed = TRUE),function(x) x[1]))
  head(tsne_res)
  
  ## 拟杆菌Bacteroides、普氏菌Prevotella和瘤胃球菌Ruminococcus
  
  ggplot(tsne_res,aes(tSNE1,tSNE2,color=tsne_res$group)) + 
    geom_point() + theme_bw() + 
    #scale_color_gradient(low = "blue",high = "red",) +
    #geom_hline(yintercept = 0,lty=2,col="red") + 
    #geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
    #stat_ellipse(aes(tSNE1,tSNE2, color = tsne_res$group), geom = "polygon", alpha = 0.3, levels = 0.99)+
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(title = "MT-RNA")
  
}

#summarize GO level matrix
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/")
  #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/GSE136243_gene_TPM_all_samples.txt",sep="\t",header = T, row.names = 1)
  #counts <- read.csv("/Users/yuhuan/Desktop/PBMC/GSE120663_TPM.csv",header = TRUE, row.names = 1)
  #counts <- counts[,-1] # for GSE40174, whose 1st colum is ENTREZGENE_ID
  #counts <- counts[,-c(1,2)] # for GSE136243, whose 1st colum is gene_name, 2nd colum is gene_biotype
  
  counts <- read.csv("/Users/yuhuan/Desktop/PBMC/pico_PBMC_TPM.csv",header = TRUE, row.names = 1)
  
  counts <- log2(as.matrix(counts)+1)
  
  # set progress
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = length(unique(GO_ID_NAME$TERM)), clear = FALSE, width= 60) 
  #get all genes in each pathway in KEGG
  n=1
  pathway_count={}
  while(n <= length(unique(GO_ID_NAME$TERM))){
    pathway_name <- as.character(unique(GO_ID_NAME$TERM)[n])
    pathway <- subset(GO_ID_NAME,TERM==pathway_name)
    
    candidate <- pathway$ensembl_gene_id #for colnames with ensemble id
    #candidate <- pathway$hgnc_symbol #for colnames with gene symbol
    
    candidate <- unique(candidate)
    
    #summary one pathway total count in each sample and make 1 row dataframe
    j=1
    pathway_gene_count={}
    while(j<=length(candidate)){
      target <- candidate[j]
      if(length(grep(target,rownames(counts)))==0) {
        #print(paste0("No ",target," in this dataset."))
        j=j+1
      } else {
        #temp <- counts[which(rownames(counts)==target),]
        temp <- counts[grep(target,rownames(counts),fixed=TRUE),]
        pathway_gene_count <- rbind(pathway_gene_count,temp)
        j=j+1
      }
    }
    if(is.null(pathway_gene_count)){
      n=n+1
      next;
    } else {
      one_pathway <- as.data.frame(t(colMeans(pathway_gene_count))) #colSums
      rownames(one_pathway) <- pathway_name
      pathway_count <- rbind(pathway_count,one_pathway)
      n=n+1
    }
    pb$tick()
    Sys.sleep(1 / 100)
    
  }
  write.csv(pathway_count,"pico_CRC_T_N_PBMC_plasma_TPM_GO_BP_level_log2.csv")
}

#volcanoplot for pico PKU
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome")
  mat_raw <- read.table("Raw_counts_passed_QC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
  des <- read.csv("des_CancervsNC_PKU.csv", header = TRUE, check.names=FALSE, sep=',')
  
  samples <- des$samples
  #samples <- gsub(".","-",samples,fixed=TRUE)
  #samples <- gsub("X","",samples,fixed=TRUE)
  group <- des$group
  batch <- des$batch
  
  i=1
  mat=as.data.frame(array(dim=c(length(rownames(mat_raw)),1)))
  while (i<=length(samples)) {
    temp <- mat_raw[,which(colnames(mat_raw)==samples[i])]
    #temp <- mat_raw[,grep(samples[i],colnames(mat_raw),fixed=TRUE)]
    #print(i)
    #print(which(colnames(mat_raw)==samples[i]))
    temp <- as.data.frame(temp)
    colnames(temp) <- samples[i]
    mat <- cbind(mat,temp)
    i=i+1
  }
  
  mat <- mat[,-1]
  rownames(mat) <- rownames(mat_raw)
  
  y <- DGEList(counts=mat, samples=samples, group=group)     #规定输入格式，必须包含这几个参数
  y <- calcNormFactors(y, method="TMM")                      #归一化处理
  y <- estimateDisp(y)
  test <- exactTest(y, pair = c("negative","positive"), dispersion = "auto")
  #广义线性模型可以针对多因素的情况，在这个场景下，只比较cancer与normal用exactTest就足够了
  {
    #design <- model.matrix(~0+group)   
    #y <- estimateDisp(y, design)          #广义线性模型计算离散度（common&trended&tagwise）
    #fit <- glmFit(y, design)              #差异表达分析函数                     
    #test <- glmLRT(fit, coef=2)           #差异表达分析函数 
  }
  
  
  
  res <- topTags(test, n=nrow(mat), sort.by='none')     #输出计算差异基因的结果
  res <- cbind(res$table, baseMean=2^(res$table$logCPM)) #输出差异基因结果加上baseMean一栏
  mapped_names <- colnames(res)                   #后面都是整理+输出
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  write.table(res, "/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/differential expression by edgeR exact/edger_exact_pico_PKU.csv", sep=',', quote=FALSE, row.names=TRUE) #输出文件，规定名字
}
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/differential expression by edgeR exact/")
  volcano <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/differential expression by edgeR exact/edger_exact_pico_PKU.csv", header = T, row.names = 1, fileEncoding = 'utf-8', fill = T)     #导入文件
  volcano$threshold <- as.factor(ifelse(volcano$padj < 0.05 & abs(volcano$log2FoldChange) >=0,ifelse(volcano$log2FoldChange > 0 ,'Up','Down'),'Not'))   #设置标记标签(factor类型）
  write.csv(volcano,"pico_PKU_CancersvsNC.csv")    #输出文件
  head(volcano)
  #volcano <- volcano[complete.cases(volcano), ]
  #head(volcano)
  #rownames(volcano) <- volcano[,1]
  volcano_plot <- ggplot(data=volcano, aes(x=log2FoldChange, y =-log10(padj), colour=threshold)) + #颜色根据threshold分（是factor）
    scale_color_manual(values=c("blue", "grey","red"))+ #设置颜色
    geom_point() + #画散点
    #xlim(c(-6,8.5)) +
    #ylim(c(0,5)) +
    theme_bw(base_size = 12, base_family = "Arial") + #装饰
    geom_vline(xintercept=c(-0,0),lty=4,col="grey",lwd=0.6,aes(0,4))+ #阀线
    geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6,aes(-3,6))+ #阀线
    theme(#legend.position="right", #装饰（下面都是）
      legend.position="right",
      panel.grid=element_blank(),
      #legend.title = element_blank(),
      #legend.text= element_text(face="bold", color="black",family = "Times", size=8),
      plot.title = element_text(hjust = 0.5,size=30,face="bold"),
      axis.text.x = element_text(face="bold",family = "Arial",color="black", size=30),
      axis.text.y = element_text(face="bold",family = "Arial",color="black", size=30),
      axis.title.x = element_text(face="bold",family = "Arial",color="black", size=36),
      axis.title.y = element_text(face="bold",family = "Arial",color="black", size=36))+
    labs(x="log2(fold change)",y="-log10 (p.adj)",title="Different Gene Expression", face="bold")
  volcano_plot   
  label <- read.csv("pico_PKU_CancersvsNC_label.csv",sep=",",header=T)
  head(label,20)
  volcano_plot+geom_text_repel(aes(x=label$log2FoldChange, y =-log10(label$padj),
                                   label = label$Gene, vjust=0.5,hjust=0.2),size=8,color="black")
}

##KEGG
KEGG <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/differential expression by edgeR exact/KEGG_pico_PKU.txt",header=T,sep="\t",quote = "")

KEGG$Term_short <- as.character(lapply(strsplit(as.character(KEGG$Term),":",fixed = TRUE),function(x) x[2]))

KEGG <- head(KEGG[order(KEGG$PValue,decreasing = FALSE),],10)

KEGG$Term_short <- factor(KEGG$Term_short,levels=KEGG[order(KEGG$PValue,decreasing = TRUE),]$Term_short)

ggplot(KEGG,aes(x=Fold.Enrichment,y=Term_short))+
  geom_point(aes(size=Count,color=-1*log10(PValue)))+
  scale_colour_gradient(low="green",high="red")+
  labs(color=expression(-log[10](P.value)),
       size="Gene number",
       x="Fold enrichment")+
  theme_bw()+
  theme(axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
        axis.title.x = element_text(size=rel(1.3)),
        axis.title.y = element_blank())

##correlation with Age
{
  Age <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/Age_pico_PKU.csv", sep = ",", header =TRUE, row.names = 1)  #输入文件1
  
  r <- cor(Age$Age,Age$Ribosome,method="pearson")       #相关性分析（皮尔森）也就是后面标注的R相关系数
  ggplot(data=Age, aes(x=Age, y=Ribosome)) +           #后面就是画散点图分析来相关性
    geom_point(color="#d7191c") +
    #geom_smooth(method="lm",color="#1a9641") +
    geom_text(aes(x=70, y=2,label=paste("R","=",signif(r,3),seq="")),color="#fdae61",size=9)+ 
    #这一步有标识R的一步：用geom_text添加，设置包括：位置、和拼接内容，颜色，字体大小
    theme_bw()+
    xlab("Age")+
    ylab("Ribosome pathway")+
    theme(
      axis.title = element_text(face="bold", color="black", size=24),
      axis.text = element_text(face="bold",  color="black", size=24)
    )
  }

##plot PBMC probes
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/")
  counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/GSE27562.txt",sep="\t",header = T, row.names = 1)
  counts <- log2(as.matrix(counts)+1)
  RP_array <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/RP-RNA-array_U133_Plus_2.csv",header = T)
  
  # set progress
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = length(unique(RP_array$DESCRPTION)), clear = FALSE, width= 60) 
  #get all genes in each pathway in KEGG
  n=1
  pathway_count={}
  while(n <= length(unique(RP_array$DESCRPTION))){
    pathway_name <- as.character(unique(RP_array$DESCRPTION)[n])
    pathway <- subset(RP_array,DESCRPTION==pathway_name)
    
    candidate <- pathway$U133_Plus_2 #for colnames with ensemble id
    #candidate <- pathway$hgnc_symbol #for colnames with gene symbol
    
    candidate <- unique(candidate)
    
    #summary one pathway total count in each sample and make 1 row dataframe
    j=1
    pathway_gene_count={}
    while(j<=length(candidate)){
      target <- candidate[j]
      if(length(grep(target,rownames(counts)))==0) {
        #print(paste0("No ",target," in this dataset."))
        j=j+1
      } else {
        #temp <- counts[which(rownames(counts)==target),]
        temp <- counts[grep(target,rownames(counts),fixed=TRUE),]
        pathway_gene_count <- rbind(pathway_gene_count,temp)
        j=j+1
      }
    }
    if(is.null(pathway_gene_count)){
      n=n+1
      next;
    } else {
      one_pathway <- as.data.frame(t(colMeans(pathway_gene_count))) #colSums
      rownames(one_pathway) <- pathway_name
      pathway_count <- rbind(pathway_count,one_pathway)
      n=n+1
    }
    pb$tick()
    Sys.sleep(1 / 100)
    
  }
  write.csv(pathway_count,"GSE27562_RNA_TPM_pathway_level_log2.csv")
}


##deconvolution by cibersortx

{
  #install.packages("remotes")
  #remotes::install_github("icbi-lab/immunedeconv")
  library(ggplot2)
  library(immunedeconv)
  library(tidyverse)
}
##GTEx
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/GTEx/TSS/")

signature <- read.csv("signature_matrix_TSS0.9_TPM_top20.txt",sep="\t",header = TRUE,row.names = 1)
col_annotation <- read.csv("sample_ids_annotation.txt",sep = "\t", header = FALSE,row.names = 1)
row_annotation <- read.csv("signature_genes_forplot.txt", sep = "\t", header = FALSE, row.names = 1)
colnames(row_annotation) <- c("Tissue_specific_gene")
colnames(col_annotation) <- c("Tissue_type","Tissue_type_specific")
rownames(col_annotation)<- gsub("-",".",rownames(col_annotation))

#signature matrix plot
bk = unique(c(seq(-0.5,0.5, length=100)))
pheatmap(
  signature,
  breaks = bk,
  annotation_col = col_annotation,
  annotation_row = row_annotation,
  scale = "row",
  cluster_cols = FALSE,cluster_rows = FALSE,
  show_colnames=FALSE, 
  show_rownames=FALSE,
  #cluster_cols = FALSE,
  colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 5)

signature_matrix <- read.csv("signature_matrix_TSS0.9_TPM_top20_median.txt",sep = "\t", header = T,row.names = 1)
pheatmap(
  signature_matrix,
  breaks = bk,
  #annotation_col = col_annotation,
  annotation_row = row_annotation,
  scale = "row",
  cluster_cols = FALSE,cluster_rows = FALSE,
  show_colnames=TRUE, 
  show_rownames=FALSE,
  #cluster_cols = FALSE,
  colorRampPalette(c("blue","white","red"))(100),
  fontsize_col = 10)


#preprocess for cibersortx
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/GTEx/forcibersortx_TSS/")
signature <- read.csv("signature_genes_forplot.txt",sep="\t",header = FALSE,row.names = 1)
mixture <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/pico_tissue_TPM.txt",sep="\t",header = T, row.names = 1)

#rownames(mixture) <- as.character(lapply(strsplit(rownames(mixture),".",fixed = TRUE),function(x) x[1]))
rownames(signature) <- as.character(lapply(strsplit(rownames(signature),".",fixed = TRUE),function(x) x[1]))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
annotations <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), 
                     filters = "ensembl_gene_id", 
                     values=rownames(signature), mart= mart,useCache = FALSE)



#get RNAs
j=1
pathway_gene_count={}
while(j<=nrow(annotations)){
  target <- annotations[j,2]
  gene_symbol <- annotations[j,1]
  if(length(grep(target,rownames(mixture)))==0) {
    #print(paste0("No ",target," in this dataset."))
    temp <- as.data.frame(array(,dim=c(1,ncol(mixture))))
    temp[1,] <- 0
    rownames(temp) <- target
    colnames(temp) <- colnames(mixture)
  } else {
    #temp <- mixture[which(rownames(mixture)==target),]
    #temp <- mixture[grep(target,rownames(mixture),fixed=TRUE),]  #for ensg
    temp <- mixture[grep(target,rownames(mixture),fixed=TRUE),] #for gene symbol
    rownames(temp) <- target
  }
  pathway_gene_count <- rbind(pathway_gene_count,temp)
  j=j+1
}

#pathway_gene_count$symbol <- annotations$hgnc_symbol
write.csv(pathway_gene_count,"pico_tissue_forGTEx_signature_TSS.csv")
#write.csv(signature,"signature_matrix_DE_top50_TPM_rmdup_ENSG.csv",quote = F)
}

#Cibersort_locally
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/GTEx/forcibersortx/")
source('Cibersort.R')  
res_cibersort <- CIBERSORT("signature_matrix_DE_top50_TPM_rmdup_ENSG_median.txt","Plasma_forGTEx_signature.txt",perm = 1000, QN = F)
write.table(res_cibersort, file="res_cibersort.txt", sep="\t", col.names=T, row.names=F, quote=F)

library(ggplot2)
library(tidyverse)
colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

my_theme <- function(){
  theme(panel.grid = element_blank(),       # 网格线
        panel.border = element_blank(),     # 面板边框
        legend.position="right",            # legend位置
        legend.text = element_text(size=8), # legend内容大小
        legend.title = element_text(size=8),# legend标题大小
        axis.line = element_line(size=1),   # 坐标轴线
        text = element_text(family="Times"),# 文本字体
        axis.text.y = element_text(size = 8,face='bold',color='black'),# y轴标签样式
        axis.text.x = element_text(size = 8,face='bold',color='black',angle=90,hjust=1),        # x轴标签样式，angle=45 倾斜 45 度
        axis.title = element_text(size=10,face="bold"),  # 轴标题
        plot.title = element_text(hjust=0.5,size=10))    # 距，设置绘图区域距离边的据类，上、右、下、左
}  

p1 <- res_cibersort[,1:22] %>% reshape2::melt() %>%
  ggplot(aes(x=Var1,y=value,fill=Var2)) +
  geom_bar(stat='identity') +
  coord_flip()  +
  scale_fill_manual(values =colour ) +
  theme_bw()+ theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+my_theme()

pdf("cibersort.pdf")
p1
}

#stackplot
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/LM22/LM22_deconvolution/PBMC/")
  composition <- read.csv("tissue_origin_PBMC.csv",header = TRUE)
  
  stackplot <- melt(composition)
  
  #boxplot
  {
  stackplot$Type <- as.character(lapply(strsplit(as.character(stackplot$variable),".",fixed = TRUE),function(x) x[3])) 
  stackplot$Type <- gsub("T","Tumor Tissue", stackplot$Type)
  stackplot$Type <- gsub("N","Normal Tissue", stackplot$Type)
  
  my_comparisons <- list(c("Normal Tissue","Tumor Tissue"))
  p <- ggplot(stackplot[which(stackplot$Average=="Blood-Vessel"),],aes(x=Type,y=value,fill=Type))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
    geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
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
      axis.title.y = element_text(face="bold",color="black", size=24))
  
  m=-1
  if(m>0){
    p <- p+stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test",
                              method.args = list(alternative = "greater"),
                              label = "p.signif"
    )+labs(x="",y="Fraction (%)",title="wilcox.test.greater", face="bold",fill="Type")
    p
  } else if(m==0){
    p <- p+stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test",
                              method.args = list(alternative = "two.sided"),
                              label = "p.signif"
    )+labs(x="",y="Fraction (%)",title="wilcox.test.twoside", face="bold",fill="Type")
    p
  } else {
    p <- p+stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test",
                              method.args = list(alternative = "less"),
                              label = "p.signif"
    )+labs(x="",y="Fraction (%)",title="Blood-Vessel\nwilcox.test.less", face="bold",fill="Type")
    p
  }
}
  
  #plot plasma
  stackplot$Average <- factor(stackplot$Average,levels = c("Blood-Whole-Blood","Blood-Vessel","Colon","Esophagus-Muscularis","Liver","Lung","Stomach"))
  stackplot$variable <- factor(stackplot$variable,levels = c("Colorectum","Stomach","Lung","Liver","Esophagus","Healthy"))
  
  #plot tissue
  comp_tmp <- stackplot[which(stackplot$Average=="T cells CD8"),]
  
  stackplot$Average <- factor(stackplot$Average,levels = unique(sort(stackplot$Average)))
  stackplot$variable <- factor(stackplot$variable,levels = comp_tmp[order(comp_tmp$value,decreasing = TRUE),2])
  
  ggplot(stackplot,aes(x=stackplot$variable,y=stackplot$value*100,fill = stackplot$Average)) + geom_bar(stat = "identity", width=0.5, col='black') +
    theme_bw()+
    theme(#legend.position="bottom",
      legend.position="right",
      panel.grid=element_blank(),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20,angle = 90,hjust = 1,vjust =0.5),
      axis.text.y = element_text(face="bold",  color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    ylab("Fraction(%)")+
    xlab("")+
    labs(fill="Tissue")+
    #geom_vline(aes(xintercept=6.5))+
    scale_y_continuous(breaks = c(0,25,50,75,100),labels = c("0","25","50","75","100"),expand = c(0,0),limits = c(0,103))
    scale_fill_nejm(alpha = 1)
}

#make median matrix and tissue specific score
{
  sample_info <- as.data.frame(list(c("a","b","c","d","e","f","g","h","i"),c("A","A","A","B","B","B","C","C","C")))
  colnames(sample_info) <- c("sample","group")
  mat <- as.data.frame(list(c(100,10,1),c(110,8,2),c(120,6,3),c(10,110,1),c(8,120,2),c(6,100,3),c(1,10,110),c(2,8,120),c(3,6,100)))
  colnames(mat) <- c("a","b","c","d","e","f","g","h","i")
  rownames(mat) <- c("gene1","gene2","gene3")
  
  group <- unique(sort(sample_info$group))
  i=1
  while(i<=length(group)){
    sample <- sample_info[which(sample_info$group==group[i]),]$sample
    median <- as.data.frame(rowMedians(as.matrix(mat[,sample])))
    colnames(median) <- group[i]
    rownames(median) <- rownames(mat)
    if(i==1){
      median_matrix <- median
    } else {
      median_matrix <- cbind(median_matrix,median)
    }
    i=i+1
  }
  
  j=1
  sum <- as.data.frame(rowSums(median_matrix))
  while(j<=ncol(median_matrix)){
    TSS <- median_matrix[,j]/sum
    colnames(TSS) <- colnames(median_matrix)[j]
    if(j==1){
      TSS_matrix <- TSS
    } else {
      TSS_matrix <- cbind(TSS_matrix,TSS)
    }
    j=j+1
  }
  
  }

#LM22
#input mixture pre-process
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/")
  mixture <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/pico_tissue_TPM.txt",sep="\t",header = T, row.names = 1)
  signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Deconvolution/LM22.txt",sep="\t",header = T, row.names = 1)
  
  #rownames(mixture) <- as.character(lapply(strsplit(rownames(mixture),".",fixed = TRUE),function(x) x[1])) 
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  annotations <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), 
                       filters = "hgnc_symbol", 
                       values=rownames(signature), mart= mart,useCache = FALSE)
  
  #get RNAs
  j=1
  pathway_gene_count={}
  while(j<=nrow(annotations)){
    target <- annotations[j,2]
    gene_symbol <- annotations[j,1]
    if(length(grep(target,rownames(mixture)))==0) {
      #print(paste0("No ",target," in this dataset."))
      temp <- as.data.frame(array(,dim=c(1,ncol(mixture))))
      temp[1,] <- 0
      rownames(temp) <-target
      colnames(temp) <- colnames(mixture)
    } else {
      #temp <- mixture[which(rownames(mixture)==target),]
      #temp <- mixture[grep(target,rownames(mixture),fixed=TRUE),]  #for ensg
      temp <- mixture[grep(target,rownames(mixture),fixed=TRUE),] #for gene symbol
      rownames(temp) <- target
    }
      pathway_gene_count <- rbind(pathway_gene_count,temp)
      j=j+1
  }
  
  pathway_gene_count$symbol <- annotations$hgnc_symbol
  write.csv(pathway_gene_count,"Tissue_forLM22.csv")
  }



###RP-RNA features

#SNP
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/RP-RNA features/SNP")
{
#read in matirx and transpose (because rownames cannot have repeat values, by colnames can)
SNP <- read.csv("SNP_matrix_with_annotation.txt",sep = "\t",header = TRUE)
SNP <- lapply(SNP, function(x) gsub("0.5","1",x))
SNP<- as.data.frame(SNP)
SNP <- t(SNP)
colnames(SNP) <- SNP[1,]
SNP <- SNP[-1,]
SNP[is.na(SNP)] <- 0 
}


SNP_numeric <- sapply(as.data.frame(SNP), function(x) {as.numeric(as.character(x))})

SNP_genevote <- rowsum(t(SNP_numeric), group = colnames(SNP_numeric), na.rm = T)

colnames(SNP_genevote) <- rownames(SNP)

write.csv(SNP_genevote,"SNP_variationlevel_count.csv")
}


#merge martrix
{
plasma <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/Raw_counts_passed_QC.txt",sep = "\t",header =T)
PBMC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Ribosome/PBMC.txt",sep = "\t",header =T)
tissue <- read.csv("/Users/yuhuan/Desktop/tissue/pico_tissue_rawcount.txt",sep = "\t",header =T)
PBMC$ENSG <- as.character(lapply(strsplit(as.character(PBMC$feature),"\\|"),function(x) x[1]))
tissue$ENSG <- as.character(lapply(strsplit(as.character(tissue$Geneid.Length),"\\|"),function(x) x[1]))
plasma$ENSG <- as.character(lapply(strsplit(as.character(plasma$feature),"\\|"),function(x) x[1]))

tissue_PBMC <- full_join(tissue,PBMC,by = c("ENSG"="ENSG"))
tissue_PBMC_plasma <- full_join(tissue_PBMC,plasma,by= c("ENSG"="ENSG"))
tissue_PBMC_plasma[is.na(tissue_PBMC_plasma)] <- 0
colnames(tissue_PBMC_plasma) <- gsub(".x",".PBMC",colnames(tissue_PBMC_plasma))
colnames(tissue_PBMC_plasma) <- gsub(".y","",colnames(tissue_PBMC_plasma))
write.table(tissue_PBMC_plasma,"/Users/yuhuan/Desktop/tissue/tissue_PBMC_plamsa.txt",sep = "\t",quote = FALSE)
}




