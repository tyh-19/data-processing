#library preparation
{
  library(reshape2)
  library(edgeR) #for RNA expression and DNA methylation differential analysis
  library(ggplot2) #for plot
  library(sunburstR) #sunburst plot
  library(htmlwidgets) #for html sunburst plot
  library(extrafont) #for plot font load
  fonts() #load fonts
  library(ggsci) #for plot color theme 
  library(cowplot) #for plot arrange
  library(tidyverse) #for dataframe process
  library(progress) #for progress bar
  library(dplyr) # for data manipulation
  library(clusterProfiler) #for enrichment analysis
  library(biomaRt) #for gene annotations
  library(multiMiR) #for miRNA target gene predict
  library(ASEP) #for allele specific expression differential analysis at gene level
  library(ggpubr) #for ggarrange
  library(pspearman) #for spearman correlation
  library(glmnet) #for cv.glmnet
  library(randomForest) #for tuneRF
}

#function
{
  #sigmoid
  sigmoid = function(x) {
    1 / (1 + exp(-x))
  }
  
  ##preparation
  library(pspearman)
  spearman_CI <- function(x, y, alpha = 0.05){
    rs <- cor(x, y, method = "spearman", use = "complete.obs")
    n <- sum(complete.cases(x, y))
    sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
  }
  
  #enmuerate correlation calculation
  enmuerate_correlation <- function(forcor,cor_methods,output_dir,max = ncol(forcor)-1){
    ##initiation
    #output
    #output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation/"
    #method
    #cor_methods <- "spearman"
    #read in data matrix, row: sample, column: numeric score. 
    #forcor <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation/Correlation between plasma expression and PBMC immune fraction.csv",header = TRUE, row.names = 1)
    #forcor <- forCorrlation
    dir.create(output_dir)
    message("Caculate R for: ",max,"(Number of interested subjects)")
    
    pb <- progress_bar$new(
      format = "  Processing [:bar] :percent eta: :eta",
      total = max, clear = FALSE, width= 60)
    
    result_final <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result_final) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    
    result <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    i=1
    while(i<= max){
      message(colnames(forcor[,i,drop = FALSE]))
      j=1
      while(j<=(ncol(forcor)-i)){
        message(colnames(forcor[,i,drop = FALSE])," vs ",colnames(forcor[,j+i,drop = FALSE]))
        test <- data.frame("A"=forcor[,i],"B"=forcor[,j+i])
        test[test=="x"] <- NA
        test <- na.omit(test)
        if(nrow(test)==0) {
          result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),"Not available","Not available","Not available","Not available","Not available","Not available")
          colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
          result <- rbind(result,result_tmp)
          j=j+1
        } else {
          if(cor_methods=="pearson"){
            r_twosided <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods)
            #r_twosided <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods)
            if (is.na(r_twosided$estimate)) {
              r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "two.sided")
              #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
            } else if (r_twosided$estimate>0){
              r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "greater")
              #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "greater")
            } else if (r_twosided$estimate<0) {
              r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "less")
              #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "less")
            } else {
              r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "two.sided")
              #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
            }
            result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,r$conf.int[1],r$conf.int[2],r_twosided$p.value)
          } else if(cor_methods=="spearman") {
            r_twosided <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
            #r_twosided <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
            if (is.na(r_twosided$estimate)) {
              r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
              #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
            } else if(r_twosided$estimate>0){
              r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "greater", approximation = "exact")
              #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "greater", approximation = "exact")
            } else if (r_twosided$estimate<0) {
              r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "less", approximation = "exact")
              #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "less", approximation = "exact")
            } else {
              r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
              #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
            }
            result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,spearman_CI(as.numeric(test$A),as.numeric(test$B))[1],spearman_CI(as.numeric(test$A),as.numeric(test$B))[2],r_twosided$p.value)
            #result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,spearman_CI(forcor[,i],forcor[,i+j])[1],spearman_CI(forcor[,i],forcor[,i+j])[2],r_twosided$p.value)
          }
          
          colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
          result <- rbind(result,result_tmp)
          j=j+1
        }
      }
      
      if(i%%20==0){
        result_final <- rbind(result_final,result)
        write.csv(result,paste0(output_dir,i/20,"_Result20.csv"))
        result <- as.data.frame(matrix(numeric(0),ncol=7))
        colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
        i=i+1
      } else if(i==(ncol(forcor)-1)){
        result_final <- rbind(result_final,result)
        write.csv(result,paste0(output_dir,i%/%20+1,"_Result",i%%20,".csv"))
        i=i+1
      } else {
        i=i+1
      }
      pb$tick()
      Sys.sleep(1 / 100)
    }
    write.csv(result_final,paste0(output_dir,"Result_final.csv"))
    
    result_final2 <- as.data.frame(matrix(numeric(0),ncol = ncol(forcor),nrow = ))
  }
  
  #paired correlation
  paired_correlation <- function(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,pvalue_cutoff=1){
    message("Please make sure input column are samples ids, line are gene ids.")
    forcor1 <- t(forcor1[gene_ids,sample_ids])
    forcor2 <- t(forcor2[gene_ids,sample_ids])
    pb <- progress_bar$new(
      format = "  Processing [:bar] :percent eta: :eta",
      total = ncol(forcor1), clear = FALSE, width= 60)
    
    result_final <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result_final) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    
    result <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    
    j=1
    while(j<=(ncol(forcor1))){
      #message(j)
      #message(colnames(forcor[,i,drop = FALSE])," vs ",colnames(forcor[,j+i,drop = FALSE]))
      test <- data.frame("A"=forcor1[,j],"B"=forcor2[,j])
      test <- test[order(test$A),]
      #test[test=="x"] <- NA
      test <- na.omit(test)
      test$ID <- factor(rownames(test),levels = rownames(test))
      
      if(nrow(test)==0) {
        result_tmp <- data.frame(paste0(colnames(forcor1)[j]," vs ",colnames(forcor2)[j]),"Not available","Not available","Not available","Not available","Not available","Not available")
        colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
        result <- rbind(result,result_tmp)
        j=j+1
      } else {
        if(cor_methods=="pearson"){
          r_twosided <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods)
          #r_twosided <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods)
          if (is.na(r_twosided$estimate)) {
            r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "two.sided")
            #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
          } else if (r_twosided$estimate>0){
            r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "greater")
            #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "greater")
          } else if (r_twosided$estimate<0) {
            r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "less")
            #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "less")
          } else {
            r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "two.sided")
            #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
          }
          result_tmp <- data.frame(paste0(colnames(forcor1)[j]," vs ",colnames(forcor2)[j]),r$estimate,r$p.value,r$method,r$conf.int[1],r$conf.int[2],r_twosided$p.value)
        } else if(cor_methods=="spearman") {
          r_twosided <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
          #r_twosided <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
          if (is.na(r_twosided$estimate)) {
            r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
            #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
          } else if(r_twosided$estimate>0){
            r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "greater", approximation = "exact")
            #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "greater", approximation = "exact")
          } else if (r_twosided$estimate<0) {
            r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "less", approximation = "exact")
            #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "less", approximation = "exact")
          } else {
            r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
            #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
          }
          result_tmp <- data.frame(paste0(colnames(forcor1)[j]," vs ",colnames(forcor2)[j]),r$estimate,r$p.value,r$method,spearman_CI(as.numeric(test$A),as.numeric(test$B))[1],spearman_CI(as.numeric(test$A),as.numeric(test$B))[2],r_twosided$p.value)
          #result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,spearman_CI(forcor[,i],forcor[,i+j])[1],spearman_CI(forcor[,i],forcor[,i+j])[2],r_twosided$p.value)
        }
        
        colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
        result <- rbind(result,result_tmp)
        
        if(j%%5000==0){
          result_final <- rbind(result_final,result)
          write.csv(result,paste0(output_dir,j/5000,"_Result5000.csv"))
          result <- as.data.frame(matrix(numeric(0),ncol=7))
          colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
          j=j+1
        } else if(j==(ncol(forcor1)-1)){
          result_final <- rbind(result_final,result)
          write.csv(result,paste0(output_dir,j%/%5000+1,"_Result",j%%5000,".csv"))
          j=j+1
        } else {
          j=j+1
        }
        
      if(is.na(r$p.value)){next}
      if(r$p.value<=pvalue_cutoff){
          p1 <-
            ggplot(test,aes(x=ID,y=0,fill= as.numeric(as.factor(test$A))))+
            geom_tile()+
            scale_fill_gradient2(low="black",mid = "white", high = "red", midpoint = length(test$A)/2)+
            #scale_fill_manual(values=c("-"="white","early"="sky blue","late"="red"))+
            theme_bw()+
            theme(#legend.position="right",
              plot.margin = unit(x=c(0,0,0,40),units="pt"),
              legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
              legend.title = element_blank(),
              #legend.text= element_blank(),
              plot.title = element_blank(),
              axis.ticks = element_blank(),
              #axis.text.x = element_blank(),
              axis.line.y = element_blank(),
              axis.text.x= element_blank(),
              #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
              axis.text.y = element_text(face="bold", color="black", size=0, angle = 90,hjust=1,vjust = 0.5),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
          
          p2 <-
            ggplot(test,aes(x=ID,y=0,fill= as.numeric(as.factor(test$B))))+
            geom_tile()+
            scale_fill_gradient2(low="black",mid = "white", high = "red", midpoint = length(test$A)/2)+
            theme_bw()+
            theme(
              plot.margin = unit(x=c(0,0,0,40),units="pt"),
              legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
              legend.title = element_blank(),
              plot.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line.y = element_blank(),
              #axis.text.x= element_blank(),
              axis.text.x = element_text( color="black", size=15, angle = 45,hjust = 1,vjust = 1),
              axis.text.y = element_text(face="bold", color="black", size=0, angle = 90,hjust=1,vjust = 0.5),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
          
          p_final <- ggarrange(p1,p2,
                               ncol = 1,nrow = 2,align = "v",heights = c(1,2.2),legend = "none")
          plot_name <- gsub("/","-",colnames(forcor1)[j-1])
          ggsave(p_final, path = output_dir, filename = paste0(r$estimate,"_",plot_name,".pdf"), device = "pdf",width = 8.14,height = 2.90)
      }
        
      }
      
      pb$tick()
      Sys.sleep(1 / 100)
    }
    write.csv(result_final,paste0(output_dir,"Paired_result_final.csv"))
  }
}

#Plot immune pathway in KEGG
{
Enrichment_STAD <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/output/Expression_STADvsNC_KEGG.txt",header = T,sep="\t",quote = "")
Enrichment_CRC <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/output/Expression_CRCvsNC_KEGG.txt",header = T,sep="\t",quote = "")
Enrichment_STAD$CancerType <- "STAD"
Enrichment_CRC$CancerType <- "CRC"
Enrichment_plot <- rbind(Enrichment_STAD,Enrichment_CRC)

Enrichment_plot <- Enrichment_plot[(Enrichment_plot$Description=="Hematopoietic cell lineage")|
                (Enrichment_plot$Description=="Complement and coagulation cascades")|
                  (Enrichment_plot$Description=="Platelet activation")|
                  (Enrichment_plot$Description=="Neutrophil extracellular trap formation")|
                  (Enrichment_plot$Description=="Toll-like receptor signaling pathway")|
                  (Enrichment_plot$Description=="Toll and Imd signaling pathway")|
                  (Enrichment_plot$Description=="NOD-like receptor signaling pathway")|
                  (Enrichment_plot$Description=="RIG-I-like receptor signaling pathway")|
                  (Enrichment_plot$Description=="Cytosolic DNA-sensing pathway")|
                  (Enrichment_plot$Description=="C-type lectin receptor signaling pathway")|
                  (Enrichment_plot$Description=="Natural killer cell mediated cytotoxicity")|
                  (Enrichment_plot$Description=="Antigen processing and presentation")|
                  (Enrichment_plot$Description=="T cell receptor signaling pathway")|
                  (Enrichment_plot$Description=="Th1 and Th2 cell differentiation")|
                  (Enrichment_plot$Description=="Th17 cell differentiation")|
                  (Enrichment_plot$Description=="IL-17 signaling pathway")|
                  (Enrichment_plot$Description=="B cell receptor signaling pathway")|
                  (Enrichment_plot$Description=="Fc epsilon RI signaling pathway")|
                  (Enrichment_plot$Description=="Fc gamma R-mediated phagocytosis")|
                  (Enrichment_plot$Description=="Leukocyte transendothelial migration")|
                  (Enrichment_plot$Description=="Intestinal immune network for IgA production")|
                  (Enrichment_plot$Description=="Chemokine signaling pathway"),]


Enrichment_plot <- Enrichment_plot[Enrichment_plot$GeneEnrichedIn=="Down regulated",]
Enrichment_plot$Label <- NA
Enrichment_plot[Enrichment_plot$p.adjust <= 0.05,]$Label <- "*"
Description_order <- aggregate(p.adjust ~ Description, Enrichment_plot, min)
Enrichment_plot$Description <- factor(Enrichment_plot$Description,levels=Description_order[order(Description_order$p.adjust,decreasing = TRUE),]$Description)
ggplot(data=Enrichment_plot)+
  geom_bar(aes(x=Description,y=Count, fill=-log10(p.adjust)),
           stat='identity') + 
  geom_text(aes(x=Description,y=Count+7,label=Label),size = 16,angle = 90)+
  facet_wrap(~Enrichment_plot$CancerType, nrow = 1)+
  coord_flip() +
  #scale_x_reverse() +
  scale_fill_gradient(expression(-log["10"](FDR)),low="dark blue", high = "red")+
  theme_light()+
  theme(panel.border = element_rect(linetype="solid",,size=0.7,colour = "black"),
        #axis.title.x=element_text(face="bold",size=12),
        #axis.text.x=element_text(face="bold",size=12),
        #axis.text.y=element_text(face="bold",size=12),
        axis.ticks.length = unit(0, "cm"),
        #axis.ticks=element_blank(),
        #axis.line = element_line(colour="black",size=1.2,linetype="solid"),
        #legend.position = "bottom",
        legend.text = element_text(size=24),
        legend.direction = NULL,
        legend.title = element_text(face="bold",size=24),
        #legend.key = element_rect(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill="white",color="transparent"),
        strip.text.x = element_text(size = 24, color = "black", face = "bold"),
        title = element_text(face="bold",size=24),
        plot.title = element_text(face="bold",size=24),
        plot.subtitle = element_text(face="bold",size=24),
        #text = element_text(size=24),
        axis.text = element_text(size=20,colour = "black"))+
  xlab("") + ylab("Gene count")
}

#siganture score boxplot
{
#make candidate list
{
  gene_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/gene/T cell receptor signaling pathway.csv",header = TRUE)
  values <- as.character(lapply(strsplit(as.character(gene_list$Gene.ID),".",fixed = TRUE),function(x) x[1]))
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  candidate_list_raw <- getBM(attributes=c("hgnc_symbol","description","ensembl_gene_id","chromosome_name","start_position","end_position","strand"), 
                              filters = "ensembl_gene_id", 
                              values=unique(values), mart= mart,useCache = FALSE)
  candidate_list_raw$strand <- gsub("-1","-",fixed = TRUE,candidate_list_raw$strand)
  candidate_list_raw$strand <- gsub("1","+",fixed = TRUE,candidate_list_raw$strand)
  candidate_list_raw <- candidate_list_raw[-grep("_",fixed=TRUE,candidate_list_raw$chromosome_name),]
  candidate_list_raw$chromosome_name <- paste0("chr",candidate_list_raw$chromosome_name)
  colnames(candidate_list_raw) <- gsub("hgnc_symbol","Gene",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("description","Description",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("ensembl_gene_id","ensembl",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("chromosome_name","chromosome",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("start_position","start",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("end_position","end",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("strand","strand",fixed = TRUE,colnames(candidate_list_raw))
  write.csv(candidate_list_raw,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/T cell receptor signaling pathway.csv",row.names = FALSE)
}
{
  gene_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/gene/GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP.csv",header = TRUE)
  values <- as.character(lapply(strsplit(as.character(gene_list$GeneName),".",fixed = TRUE),function(x) x[1]))
  gene_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/STAD_colon_coupregulated_genes.csv",header = TRUE)
  values <- as.character(lapply(strsplit(as.character(gene_list$hgnc_symbol),".",fixed = TRUE),function(x) x[1])) 
  #gene_list <- gene_list[which(gene_list$Super.Category=="Antigen presentation"),]
  
  
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  candidate_list_raw <- getBM(attributes=c("hgnc_symbol","description","ensembl_gene_id","chromosome_name","start_position","end_position","strand"), 
                              filters = "hgnc_symbol", 
                              values=unique(values), mart= mart,useCache = FALSE)
  candidate_list_raw$strand <- gsub("-1","-",fixed = TRUE,candidate_list_raw$strand)
  candidate_list_raw$strand <- gsub("1","+",fixed = TRUE,candidate_list_raw$strand)
  candidate_list_raw <- candidate_list_raw[-grep("_",fixed=TRUE,candidate_list_raw$chromosome_name),]
  candidate_list_raw$chromosome_name <- paste0("chr",candidate_list_raw$chromosome_name)
  colnames(candidate_list_raw) <- gsub("hgnc_symbol","Gene",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("description","Description",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("ensembl_gene_id","ensembl",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("chromosome_name","chromosome",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("start_position","start",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("end_position","end",fixed = TRUE,colnames(candidate_list_raw))
  colnames(candidate_list_raw) <- gsub("strand","strand",fixed = TRUE,colnames(candidate_list_raw))
  write.csv(candidate_list_raw,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/gene/STAD_colon_coupregulated_genes.csv",row.names = FALSE)
}
#Plot siganture score of immune pathways
{
pathway <- "GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP"
pathway <- "GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP"
pathway <- "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP"
pathway <- "T cell receptor signaling pathway"
pathway <- "B cell receptor signaling pathway"
pathway <- "Inhibotory immune genes"
pathway <- "PELO"
pathway <- "STAD_colon_coupregulated_genes"
i=1

prefix <- "Multiomics_Plasma"
counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/GSE174302_intron-spanning-available-samples_TPM.txt",sep = "\t", header = TRUE,check.names = FALSE)
counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/GSE133684_count_matrix_TPM.txt",sep = "\t", header = TRUE,check.names = FALSE)
counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/multiomics_plasma_TPM_intron_spanning_20220427.txt",sep = "\t", header = TRUE,check.names = FALSE, row.names = 1)
counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Plasma_TPM.txt",sep = "\t", header = TRUE,check.names = FALSE, row.names = 1)
rownames(counts) <- as.character(lapply(strsplit(as.character(rownames(counts)),".",fixed = TRUE),function(x) x[1]))
counts$Gene <- as.character(lapply(strsplit(as.character(rownames(counts)),".",fixed = TRUE),function(x) x[1]))
counts <- aggregate(counts[,-which(colnames(counts)=="Gene")], list(Gene=counts$Gene), FUN = sum)
rownames(counts) <- counts$Gene
counts <- counts[,-which(colnames(counts)=="Gene")]
#counts <- counts[,-which(colnames(counts)=="CRC-PKU-29-PBMC")]

prefix <- "Multiomics_20211113"
#counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/Expression/matrix/20211113_TPM.txt",sep = "\t", header = TRUE,check.names = FALSE)
#counts$`Gene|Length` <- as.character(lapply(strsplit(as.character(counts$`Gene|Length`),".",fixed = TRUE),function(x) x[1]))
counts <- aggregate(counts[,-1], list(Gene=counts[,1]), FUN = sum)
rownames(counts) <- counts$Gene
counts <- counts[,-1]

pathway_gene <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/gene/",pathway,".csv"),header = T)
}
#get pathway genes
{
  j=1
  pathway_gene_count={}
  while(j<=nrow(pathway_gene)){
    target <- pathway_gene[j,which(colnames(pathway_gene)=="ensembl")]
    gene_symbol <- pathway_gene[j,which(colnames(pathway_gene)=="Gene")]
    #target <- pathway_gene[j,which(colnames(pathway_gene)=="Gene.ID")]
    #gene_symbol <- pathway_gene[j,which(colnames(pathway_gene)=="Gene.Name")]
    if(length(grep(target,rownames(counts)))==0) {
      print(paste0("No ",target," in this dataset."))
      temp <- as.data.frame(array(,dim=c(1,ncol(counts))))
      temp[1,] <- 0
      rownames(temp) <- gene_symbol
      j=j+1
    } else {
      #temp <- counts[which(rownames(counts)==target),]
      temp <- counts[grep(target,rownames(counts),fixed=TRUE),]  #for ensg
      rownames(temp) <- gene_symbol
      pathway_gene_count <- rbind(pathway_gene_count,temp)
      j=j+1
    }
  }
  pathway_gene_count_log2 <- log2(pathway_gene_count+1)
  pathway_gene_count_log2_colmean <- colMeans(pathway_gene_count_log2)
  pathway_gene_count_log2_colmean <- as.data.frame(pathway_gene_count_log2_colmean)
  colnames(pathway_gene_count_log2_colmean) <- pathway
  
  pathway_gene_count_log2_colsum <- colSums(pathway_gene_count_log2)
  pathway_gene_count_log2_colsum <- as.data.frame(pathway_gene_count_log2_colsum)
  colnames(pathway_gene_count_log2_colsum) <- pathway
  
  pathway_gene_count_colmean <- colMeans(pathway_gene_count)
  pathway_gene_count_colmean <- as.data.frame(pathway_gene_count_colmean)
  colnames(pathway_gene_count_colmean) <- pathway
  
  dir.create(path = paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i]), recursive = TRUE)
  write.table(pathway_gene_count,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,"_matrix.txt"),quote = FALSE,sep = "\t")
  write.table(pathway_gene_count_log2_colmean,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,"_log2.txt"),quote = FALSE,sep = "\t")
  write.table(pathway_gene_count_log2_colsum,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,"_log2sum.txt"),quote = FALSE,sep = "\t")
  write.table(pathway_gene_count_colmean,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,".txt"),quote = FALSE,sep = "\t")
}
#plot only compares group
{
  des <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/group_info/",prefix[i],".csv"),header = TRUE, row.names = 1)
  log2_mean <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,"_log2sum.txt"),sep = "\t",header = TRUE, row.names = 1,check.names=FALSE)
  rownames(log2_mean) <- gsub(".","-",fixed= TRUE,rownames(log2_mean))
  des$sample_id <- rownames(des)
  des <- des[-which(des$group=="PBMC"),]
  #des <- des[-which(rownames(des)=="CRC-PKU-29-PBMC"),]
  log2_mean$sample_id <- rownames(log2_mean)
  boxplot <- left_join(des,log2_mean,by = c('sample_id'='sample_id'))
  
  #boxplot <- na.omit(boxplot)
  boxplot$group <- factor(boxplot$group,levels=c("HD","CRC","STAD"))
  my_comparisons <- list(c("STAD","HD"),c("CRC","HD"))
  #my_comparisons <- list(c("CRC","HD"))
  my_comparisons <- list(c("negative","positive"))
  p <- ggplot(boxplot,aes(x=group,y=boxplot[[pathway]],fill=group))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
    geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
    #scale_fill_manual(values = c("CRC"="#FCB514","positive"="red","negative"="blue")) +
    scale_fill_manual(values = c("CRC"="#FCB514","STAD"="red","HD"="blue")) +
    #scale_fill_manual(values = c("#EE7621","blue"))+
    ylab(pathway)+
    theme_bw()+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(size=1, colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=30,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
      axis.text.y = element_text(face="bold",  color="black", size=24),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
  
  dir.create(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"), recursive = TRUE)
  m=mean(boxplot[boxplot$group=="negative",ncol(boxplot)])-mean(boxplot[boxplot$group=="positive",ncol(boxplot)])
  if(m>0){
    p <- p+stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test",
                              size = 8,
                              vjust = 0.5,
                              method.args = list(alternative = "greater"),
                              label = "p.signif"
    )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.greater"), face="bold",fill="Type")
    ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_greater.pdf"),width = 6,height = 8, units = "in", dpi = 300)
  } else if(m==0){
    p <- p+stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test",
                              size = 8,
                              vjust = 0.5,
                              method.args = list(alternative = "two.sided"),
                              label = "p.signif"
    )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.twosided"), face="bold",fill="Type")
    ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_twosided.pdf"),width = 6,height = 8, units = "in", dpi = 300)
  } else {
    p <- p+stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test",
                              size = 8,
                              vjust = 0.5,
                              method.args = list(alternative = "less"),
                              label = "p.signif"
    )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.less"), face="bold",fill="Type")
    ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_less.pdf"),width = 6,height = 8, units = "in", dpi = 300)
  }
}
  
#plot compares stage in cancer for GSE174302
{
    des <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/group_info/GSE174302_CRCSTAD_vs_NC_stage.csv"),header = TRUE, row.names = 1)
    log2_mean <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,"_log2sum.txt"),sep = "\t",header = TRUE, row.names = 1,check.names=FALSE)
    rownames(log2_mean) <- gsub(".","-",fixed= TRUE,rownames(log2_mean))
    des$sample_id <- rownames(des)
    log2_mean$sample_id <- rownames(log2_mean)
    boxplot <- left_join(des,log2_mean,by = c('sample_id'='sample_id'))
    
    #boxplot <- na.omit(boxplot)
    boxplot$group <- factor(boxplot$group,levels=c("Cancer","HD"))
    #my_comparisons <- list(c("STAD","HD"),c("CRC","HD"))
    boxplot$Stage.simplified <- factor(boxplot$Stage.simplified,levels=c("No stage","Stage I","Stage II","Stage III","Stage IV"))
    my_comparisons <- list(c("Stage I","No stage"),c("Stage II","No stage"),c("Stage III","No stage"),c("Stage IV","No stage"))
    p <- ggplot(boxplot,aes(x=Stage.simplified,y=boxplot[[pathway]],fill=group))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
      scale_fill_manual(values = c("Cancer"="red","HD"="blue")) +
      #scale_fill_manual(values = c("CRC"="#FCB514","STAD"="red","HD"="blue")) +
      #scale_fill_manual(values = c("#EE7621","blue"))+
      ylab(pathway)+
      theme_bw()+
      theme(#legend.position="right",
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=30,face="bold"),
        axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))
    
    dir.create(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"), recursive = TRUE)
    m=mean(boxplot[boxplot$group=="Cancer",ncol(boxplot)])-mean(boxplot[boxplot$group=="HD",ncol(boxplot)])
    if(m>0){
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                size = 8,
                                vjust = 0.5,
                                method.args = list(alternative = "greater"),
                                label = "p.signif"
      )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.greater"), face="bold",fill="Type")
      ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_greater.pdf"),width = 6,height = 8, units = "in", dpi = 300)
    } else if(m==0){
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                size = 8,
                                vjust = 0.5,
                                method.args = list(alternative = "two.sided"),
                                label = "p.signif"
      )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.twosided"), face="bold",fill="Type")
      ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_twosided.pdf"),width = 6,height = 8, units = "in", dpi = 300)
    } else {
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                size = 8,
                                vjust = 0.5,
                                method.args = list(alternative = "less"),
                                label = "p.signif"
      )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.less"), face="bold",fill="Type")
      ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_less.pdf"),width = 9,height = 8, units = "in", dpi = 300)
    }
}
  
#plot GSE27562 (PBMC RNA array)
  {
    pathway <- "Inhibitory immune genes"
    prefix <- "GSE27562"
    pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/GSE27562_T cell receptor_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
    pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/GSE27562_B cell receptor_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
    pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/GSE27562_Inhibotory immune genes_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
    pathway_count <- as.data.frame(t(pathway_count))
    pathway_count$ID <- rownames(pathway_count)
    
    group <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/GSE27562_des.csv",header = TRUE)
    boxplot <- left_join(pathway_count,group, by=c("ID"="ID"))
    
    my_comparisons <- list(c("Cancer","Healthy"))
    boxplot$group <- factor(boxplot$group,levels = c("Healthy","Cancer"))
    p <- ggplot(boxplot,aes(x=group,y=boxplot[[pathway]],fill=group))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
      scale_fill_manual(values = c("Cancer"="red","Healthy"="blue")) +
      #scale_fill_manual(values = c("CRC"="#FCB514","STAD"="red","HD"="blue")) +
      #scale_fill_manual(values = c("#EE7621","blue"))+
      ylab(pathway)+
      theme_bw()+
      theme(#legend.position="right",
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=30,face="bold"),
        axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))
    
    dir.create(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"), recursive = TRUE)
    m=mean(boxplot[boxplot$group=="Cancer",1])-mean(boxplot[boxplot$group=="Healthy",1])
    if(m>0){
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                size = 8,
                                vjust = 0.5,
                                method.args = list(alternative = "greater"),
                                label = "p.signif"
      )+labs(x="",y="log2(Signal Intensity+1)",title=paste0(prefix[i],"\nwilcox.test.greater"), face="bold",fill="Type")
      ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_greater.pdf"),width = 6,height = 8, units = "in", dpi = 300)
    } else if(m==0){
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                size = 8,
                                vjust = 0.5,
                                method.args = list(alternative = "two.sided"),
                                label = "p.signif"
      )+labs(x="",y="log2(Signal Intensity+1)",title=paste0(prefix[i],"\nwilcox.test.twosided"), face="bold",fill="Type")
      ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_twosided.pdf"),width = 6,height = 8, units = "in", dpi = 300)
    } else {
      p <- p+stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox.test",
                                size = 8,
                                vjust = 0.5,
                                method.args = list(alternative = "less"),
                                label = "p.signif"
      )+labs(x="",y="log2(Signal Intensity+1)",title=paste0(prefix[i],"\nwilcox.test.less"), face="bold",fill="Type")
      ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_less.pdf"),width = 9,height = 8, units = "in", dpi = 300)
}
}
}
#lncRNA-GC1 and immune gene boxplot
{
  #test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/Expression/matrix/20211113_TPM.txt",sep = "\t", header = TRUE, row.names = 1)
  #plot <- as.data.frame(t(test[which(rownames(test)=="lncRNA-GC1|2145"),]))
  #test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20210722_multiomics/Altpromoter/Altpromoter_ML.txt",sep = "\t", header = TRUE, row.names = 1)
  
  #test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/pico_PBMC_featurecount_intron_spanning_noMTRNA_TPM.txt",sep = "\t", header = TRUE, row.names = 1)  
  #test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/pico_tissue_featurecount_intron_spanning_noMTRNA_TPM.txt",sep = "\t", header = TRUE, row.names = 1)  
  test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Plasma_TPM.txt",sep = "\t", header = TRUE, row.names = 1)
  
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
  
  #CD8 T cell marker
  #ENSG00000211789: TRAV12-2
  #ENSG00000153563: CD8A
  #ENSG00000172116: CD8B
  
  
  
  #CD4 memory naive T cell marker
  #"ENSG00000231160"="KLF3-AS1"
  #"ENSG00000176293"="ZNF135"
  #"ENSG00000125864"="BFSP1"
  #"ENSG00000166573"="GALR1"
  #"ENSG00000168067"="MAP4K2"
  #"ENSG00000134765"="DSC1"
  #"ENSG00000204789"="ZNF204P"
  #"ENSG00000146904"="EPHA1"
  #"ENSG00000090554"="FLT3LG"
  #"ENSG00000138795"="LEF1"
  #"ENSG00000197093"="GAL3ST4"
  #"ENSG00000154764"="WNT7A"
  #"ENSG00000129158"="SERGEF"
  #"ENSG00000164512"="ANKRD55"
  #"ENSG00000142102"="PGGHG"
  #"ENSG00000083812"="ZNF324"
  
  #CD 4 memory resting cell marker
  #"ENSG00000163346"="RCAN3"
  #"ENSG00000152518"="ZFP36L2"
  #"ENSG00000159023"="EPB41"
  #"ENSG00000134954"="ETS1"
  #"ENSG00000135722"="FBXL8"
  #"ENSG00000165496"="RPL10L"
  #"ENSG00000197635"="DPP4"
  #"ENSG00000117602"="RCAN3"
  #"ENSG00000139193"="CD27"
  #"ENSG00000170128"="GPR25"
  #"ENSG00000177272"="KCNA3"
  #"ENSG00000100100"="PIK3IP1"
  #"ENSG00000198851"="CD3E"
  #"ENSG00000107742"="SPOCK2"
  
  #gama delta
  "ENSG00000109943"="CRTAM"
  "ENSG00000158050"="DUSP2"
  "ENSG00000113088"="GZMK"
  "ENSG00000139187"="KLRG1"
  "ENSG00000147231"="RADX"
  "ENSG00000233402"="TARDBPP1"
  "ENSG00000160791"="CCR5"
  "ENSG00000198342"="ZNF442"
  "ENSG00000211829"="TRDC"
  "ENSG00000211804"="TRDV1"
  "ENSG00000211804"="TRDV2"
  
  List <-  c("ENSG00000120217"="CD274(PD-L1)",
             "ENSG00000089692"="LAG3",
             "ENSG00000109943"="CRTAM",
             "ENSG00000158050"="DUSP2",
             "ENSG00000113088"="GZMK",
             "ENSG00000139187"="KLRG1",
             "ENSG00000147231"="RADX",
             "ENSG00000233402"="TARDBPP1",
             "ENSG00000160791"="CCR5",
             "ENSG00000198342"="ZNF442",
             "ENSG00000211829"="TRDC",
             "ENSG00000211804"="TRDV1",
             "ENSG00000211804"="TRDV2",
             "ENSG00000211789"="TRAV12-2",
             "ENSG00000153563"="CD8A",
             "ENSG00000172116"="CD8B",
             "ENSG00000010610"="CD4",
             "ENSG00000178562"="CD28",
             "ENSG00000177455"="CD19",
             "ENSG00000115085"="ZAP70",
             "ENSG00000113263"="ITK",
             "ENSG00000111537"="IFNG",
             "ENSG00000074966"="RLK",
             "ENSG00000170345"="FOS",
             "ENSG00000103855"="CD276(B7-H3)",
             "ENSG00000134258"="VTCN1(B7-H4)",
             "ENSG00000114455"="HHLA2(B7-H5)",
             "ENSG00000188389"="PDCD1(PD1)",
             "ENSG00000163599"="CTLA4",
             "ENSG00000089692"="LAG3",
             "ENSG00000135077"="TIM3",
             "ENSG00000079385"="CEACAM1(TIM3 ligand)",
             "ENSG00000181847"="TIGIT",
             "ENSG00000198793"="mTOR")
  
  #List <- as.list(rownames(LM22))
  #names(List) <- rownames(LM22)
  
  i=1
  while(i <= length(List)){
    Gene_ID <- names(List)[i]
    Gene_name <- as.character(List[i])
    {
      plot <- as.data.frame(t(test[grep(Gene_ID,rownames(test)),]))
      plot$sample <- rownames(plot)
      plot <- plot[grep("pico",rownames(plot)),]
      plot <- plot[-grep("NC.PKU.mix..pico|CRC.PKU.mix1.pico|CRC.PKU.5.pico|NC.PKU.mix17.pico|STAD.PKU.4.pico",rownames(plot)),]
      
      plot$group <- as.character(lapply(strsplit(rownames(plot),".",fixed=TRUE),function(x) x[1]))
      #plot$group <- as.character(lapply(strsplit(rownames(plot),".",fixed=TRUE),function(x) tail(x,n=1)))
      
      plot$group <- gsub("NC","HD",plot$group)
      
      my_comparisons <- list(c("STAD","HD"),c("CRC","HD"))
      #my_comparisons <- list(c("CRC","HD"))
      #my_comparisons <- list(c("T","N"))
      plot$group <- factor(plot$group,levels=c("CRC","STAD","HD"))
      
      #plot$group <- gsub("CRC","GIC",plot$group)
      #plot$group <- gsub("STAD","GIC",plot$group)
      #my_comparisons <- list(c("HD","GIC"))
      #plot$group <- factor(plot$group,levels=c("HD","GIC"))
      forplot <- plot[,c(1,ncol(plot))]
      colnames(forplot) <- c("value","group")
      #p <- ggplot(forplot[-which(rownames(forplot)=="CRC.PKU.29.pico"),],aes(x=group,y=value,fill = group))+
      p <- ggplot(forplot,aes(x=group,y=value,fill = group))+
        geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
        geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
        scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue"))+
        #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
        #scale_fill_brewer(palette="Blues") +
        #ylim(0,25)+
        theme_bw()+
        xlab("")+
        ylab(Gene_name)+
        #ylab(colnames(plot)[1])+
        theme(#legend.position="right",
          legend.position="none",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(size=1, colour = "black"),
          legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
          plot.title = element_text(hjust = 0.5,size=36,face="bold"),
          axis.text.x = element_blank(),
          #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
          axis.text.y = element_text(face="bold",  color="black", size=24),
          axis.title.x = element_text(face="bold", color="black", size=24),
          axis.title.y = element_text(face="bold",color="black", size=24))+
        stat_compare_means(comparisons = my_comparisons,
                           method = "wilcox.test",
                           method.args = list(alternative = "two.sided",paired = TRUE),
                           label = "p.signif",
                           size = 10,
                           vjust = 0.5)
      ggsave(plot = p,path = "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/markers/", filename = paste0(List[i],".pdf"),device = "pdf",height = 5.0,width = 3.8)
      i=i+1
    }
  }
}

#LM22 radar plot
{
  {
    composition <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/res_cibersort_multiomics_plasma_PBMC_tissue_LM22_abs_20220608_perm1000.txt",sep="\t",header=TRUE,row.names=1)
    composition <- composition[,-which(colnames(composition)=="RMSE")]
    composition <- composition[,-which(colnames(composition)=="Correlation")]
    composition <- composition[,-which(colnames(composition)=="P.value")]
    composition <- composition[,-which(colnames(composition)=="Absolute.score..sig.score.")]
    
    composition$B.cells <- composition$B.cells.naive+composition$B.cells.memory
    composition <- composition[,-which(colnames(composition)=="B.cells.naive")]
    composition <- composition[,-which(colnames(composition)=="B.cells.memory")]
    composition$T.cells.CD4 <- composition$T.cells.CD4.naive+composition$T.cells.CD4.memory.resting+composition$T.cells.CD4.memory.activated+composition$T.cells.follicular.helper+composition$T.cells.regulatory..Tregs.
    composition <- composition[,-which(colnames(composition)=="T.cells.CD4.naive")]
    composition <- composition[,-which(colnames(composition)=="T.cells.CD4.memory.resting")]
    composition <- composition[,-which(colnames(composition)=="T.cells.CD4.memory.activated")]
    composition <- composition[,-which(colnames(composition)=="T.cells.follicular.helper")]
    composition <- composition[,-which(colnames(composition)=="T.cells.regulatory..Tregs.")]
    composition$NK.cells <- composition$NK.cells.resting+composition$NK.cells.activated
    composition <- composition[,-which(colnames(composition)=="NK.cells.resting")]
    composition <- composition[,-which(colnames(composition)=="NK.cells.activated")]
    composition$Macrophages <- composition$Macrophages.M0+composition$Macrophages.M1+composition$Macrophages.M2
    composition <- composition[,-which(colnames(composition)=="Macrophages.M0")]
    composition <- composition[,-which(colnames(composition)=="Macrophages.M1")]
    composition <- composition[,-which(colnames(composition)=="Macrophages.M2")]
    composition$Dendritic.cells <- composition$Dendritic.cells.resting+composition$Dendritic.cells.activated
    composition <- composition[,-which(colnames(composition)=="Dendritic.cells.resting")]
    composition <- composition[,-which(colnames(composition)=="Dendritic.cells.activated")]
    composition$Mast.cells <- composition$Mast.cells.resting+composition$Mast.cells.activated
    composition <- composition[,-which(colnames(composition)=="Mast.cells.resting")]
    composition <- composition[,-which(colnames(composition)=="Mast.cells.activated")]
  }
  
  #composition <- as.data.frame(t(apply(composition, 1, function(x) x / sum(as.numeric(x)) * 10^2)))
  composition <- na.omit(composition)
  composition<-composition[grep("pico",rownames(composition)),]
  composition<-composition[-grep("mix..pico",rownames(composition)),]
  composition<-composition[-grep("CRC-PKU-5-pico",rownames(composition)),]
  composition<-composition[-grep("NC-PKU-mix17-pico",rownames(composition)),]
  composition<-composition[-grep("STAD-PKU-4-pico",rownames(composition)),]
  #paired tissue and plasma
  composition<-composition[grep("NC|STAD|CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(composition)),]
  
  #composition<-composition[grep("PBMC",rownames(composition)),]
  composition<-composition[-grep("PBMC",rownames(composition)),]
  composition<-composition[-grep("pico",rownames(composition)),]
  
  
  
  composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[1]))
  composition$group<-gsub("NC","HD",composition$group)
  #composition$group<-gsub("CRC","GIC",composition$group)
  #composition$group<-gsub("STAD","GIC",composition$group)
  #composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[length(x)]))
  #composition$group<-gsub("T","Tumor",composition$group)
  #composition$group<-gsub("N","Normal",composition$group)
  
  radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=median)
  radar <- as.tibble(radar)
  x <- rev(order(radar[2,-1]))+1
  #x <- c(8,2,3,9,4,10,5,11,12,13,6,7)
  radar <- radar[,c(1,x)]
  
  colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
  colnames(radar) <- gsub("T cells CD8","CD8 T cells",colnames(radar))
  colnames(radar) <- gsub("T cells CD4","CD4 T cells",colnames(radar))
  colnames(radar) <- gsub("T cells gamma delta","γδ T cells",colnames(radar))
  colnames(radar) <- gsub("Dendritic cells","Dendritic\ncells",colnames(radar))
  #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
  ggradar(radar,grid.min = 0,grid.mid = 0.25, grid.max = 0.5,
          #axis.labels = c("B cells","Plasma cells","CD8 T cells","CD4 T cells","γδ T cells","NK cells",
          #                "Monocytes","Macrophages","Dendritic cells","Mast cells","Eosinophils","Neutrophils"),
          grid.line.width = 1,
          base.size = 58,
          values.radar = c("", "0.25", "0.5"),
          plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
          font.radar = "Arial",
          group.point.size = 3,
          group.line.width = 2,
          legend.position = "right",
          background.circle.transparency = 0.0,
          group.colours = c(CRC="#FCB514",STAD="red",HD="blue",GIC="#FCB514",Tumor="#FCB514",Normal="blue"))
}

#plot subtype
## L/R no signaificant change(may cause by too less R samples(only 5))
{
  des <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/group_info/",prefix[i],".csv"),header = TRUE, row.names = 1,check.names=FALSE)
  log2_mean <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,"_log2sum.txt"),sep = "\t",header = TRUE, row.names = 1,check.names=FALSE)
  rownames(log2_mean) <- gsub(".","-",fixed= TRUE,rownames(log2_mean))
  des$sample_id <- rownames(des)
  des <- des[-which(des$group=="PBMC"),]
  log2_mean$sample_id <- rownames(log2_mean)
  boxplot <- left_join(des,log2_mean,by = c('sample_id'='sample_id'))
  
  #boxplot <- na.omit(boxplot)
  boxplot <- boxplot[which(boxplot$group=="CRC"),]
  boxplot$`Right half / left half` <- factor(boxplot$`Right half / left half`,levels=c("L","R"))
  my_comparisons <- list(c("L","R"))
  p <- ggplot(boxplot,aes(x=boxplot$`Right half / left half`,y=boxplot$`Th1 and Th2 cell differentiation`,fill=boxplot$`Right half / left half`))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
    geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
    scale_fill_manual(values = c("#EE7621","red","blue")) +
    theme_bw()+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(size=1, colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=30,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
      axis.text.y = element_text(face="bold",  color="black", size=24),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))
  
  dir.create(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"), recursive = TRUE)
  m=0
  if(m>0){
    p <- p+stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test",
                              size = 8,
                              vjust = 0.5,
                              method.args = list(alternative = "greater"),
                              label = "p.signif"
    )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.greater"), face="bold",fill="Type")
    ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_greater.pdf"),width = 6,height = 8, units = "in", dpi = 300)
  } else if(m==0){
    p <- p+stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test",
                              size = 8,
                              vjust = 0.5,
                              method.args = list(alternative = "two.sided"),
                              label = "p.signif"
    )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.twosided"), face="bold",fill="Type")
    #ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_twosided.pdf"),width = 6,height = 8, units = "in", dpi = 300)
  } else {
    p <- p+stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test",
                              size = 8,
                              vjust = 0.5,
                              method.args = list(alternative = "less"),
                              label = "p.signif"
    )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.less"), face="bold",fill="Type")
    #ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_less.pdf"),width = 6,height = 8, units = "in", dpi = 300)
  }
  p
}

#GSEA
{
#Plot GSEA curve for immune pathways
{
  #prepare GSEA reference geneset (KEGG)
  KEGG_forGSEA <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/09.Published_data/pathway/PATH_ID_NAME_modified.csv",header  = TRUE,row.names = 1)
  GMT <- KEGG_forGSEA[,c("DESCRPTION","ensembl_gene_id")]
  
  #CRC
  {
    CRC_Differential <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/output/Expression_tissue_CRC_TvsN_edger_exact.txt",header = T,sep="\t",quote = "")
    CRC_Differential$ensemble_gene_id <- as.character(lapply(strsplit(as.character(rownames(CRC_Differential)),".",fixed = TRUE),function(x) x[1]))
    CRC_Differential <- aggregate( CRC_Differential[,-which(colnames(CRC_Differential)=="ensemble_gene_id")], list(ensemble_gene_id=CRC_Differential$ensemble_gene_id), FUN = mean)
    CRC <- CRC_Differential[grep("ENSG",CRC_Differential$ensemble_gene_id),c("ensemble_gene_id","log2FoldChange")]
    CRC_genelist <- CRC$log2FoldChange
    names(CRC_genelist) <- CRC$ensemble_gene_id
    CRC_genelist <- sort(CRC_genelist,decreasing = TRUE)
    CRC_GSEA <- GSEA(geneList = CRC_genelist,TERM2GENE = GMT, pvalueCutoff = 0.1)
    CRC_GSEA$Description
    pathway <- "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis"
    gseaplot(CRC_GSEA, geneSetID = which(CRC_GSEA$Description==pathway), title = pathway)
  }
  
  #STAD
  {
    STAD_Differential <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/output/Expression_STADvsNC_edger_exact.txt",header = T,sep="\t",quote = "")
    STAD_Differential$ensemble_gene_id <- as.character(lapply(strsplit(as.character(rownames(STAD_Differential)),".",fixed = TRUE),function(x) x[1]))
    STAD_Differential <- aggregate( STAD_Differential[,-which(colnames(STAD_Differential)=="ensemble_gene_id")], list(ensemble_gene_id=STAD_Differential$ensemble_gene_id), FUN = mean)
    STAD <- STAD_Differential[grep("ENSG",STAD_Differential$ensemble_gene_id),c("ensemble_gene_id","log2FoldChange")]
    STAD_genelist <- STAD$log2FoldChange
    names(STAD_genelist) <- STAD$ensemble_gene_id
    STAD_genelist <- sort(STAD_genelist,decreasing = TRUE)
    STAD_GSEA <- GSEA(geneList = STAD_genelist,TERM2GENE = GMT,pvalueCutoff = 0.1)
    STAD_GSEA$Description
    pathway <- "T cell receptor signaling pathway"
    gseaplot(STAD_GSEA, geneSetID = which(STAD_GSEA$Description==pathway), title = pathway)
  }
  
  #GIC
  {
    GIC_Differential <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/output/Expression_CRCvsNC_paired_withTissue_edger_exact.txt",header = T,sep="\t",quote = "")
    GIC_Differential$ensemble_gene_id <- as.character(lapply(strsplit(as.character(rownames(GIC_Differential)),".",fixed = TRUE),function(x) x[1]))
    GIC_Differential <- aggregate(GIC_Differential[,-which(colnames(GIC_Differential)=="ensemble_gene_id")], list(ensemble_gene_id=GIC_Differential$ensemble_gene_id), FUN = median)
    GIC <- GIC_Differential[grep("ENSG",GIC_Differential$ensemble_gene_id),c("ensemble_gene_id","log2FoldChange")]
    GIC_genelist <- GIC$log2FoldChange
    names(GIC_genelist) <- GIC$ensemble_gene_id
    GIC_genelist <- sort(GIC_genelist,decreasing = TRUE)
    GIC_GSEA <- GSEA(geneList = GIC_genelist,TERM2GENE = GMT, pvalueCutoff = 1, minGSSize = 20, maxGSSize = Inf,)
    GIC_GSEA$Description
    pathway <- "T cell receptor signaling pathway"
    gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
    
    pathway <- "B cell receptor signaling pathway"
    gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
    
    pathway <- "Cancer geneset (no immune or metabolism)"
    gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
    
    pathway <- "Immune, metabolism and cancer"
    gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
    
    pathway <- "Colorectal cancer"
    gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
    
    pathway <- "Gastric cancer"
    gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
    
    enrichplot::gseaplot2(GIC_GSEA,c("B cell receptor signaling pathway","T cell receptor signaling pathway","Oxidative phosphorylation"),
                          base_size = 24,
                          color = c("B cell receptor signaling pathway"="Blue",
                                    "T cell receptor signaling pathway"="dark blue",
                                    "Oxidative phosphorylation"="red"),
                          #color = "black",
                          rel_heights = c(2, 0.25, 1),
                          subplots = c(1,2,3),
                          ES_geom = "line",
                          pvalue_table = FALSE)
    library(enrichplot)
    library(gridExtra)
    library(grid)
    trace("gseaplot2", edit = TRUE)
    
    ridgeplot(GIC_GSEA)
  }
}

#Plot GSEA curve for PBMC and tissue
  {
    #prepare GSEA reference geneset (KEGG)
    KEGG_forGSEA <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/KEGG_pathway_genes/PBMC_tissue_genes的副本.csv",header  = TRUE,row.names = 1)
    GMT <- KEGG_forGSEA[,c("DESCRPTION","ensembl_gene_id")]
    
    #GIC
    {
      GIC_Differential <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/output/Expression_tissue_CRC_TvsN_edger_exact.txt",header = T,sep="\t",quote = "")
      GIC_Differential$ensemble_gene_id <- as.character(lapply(strsplit(as.character(rownames(GIC_Differential)),".",fixed = TRUE),function(x) x[1]))
      GIC_Differential <- aggregate(GIC_Differential[,-which(colnames(GIC_Differential)=="ensemble_gene_id")], list(ensemble_gene_id=GIC_Differential$ensemble_gene_id), FUN = mean)
      GIC <- GIC_Differential[grep("ENSG",GIC_Differential$ensemble_gene_id),c("ensemble_gene_id","log2FoldChange")]
      GIC_genelist <- GIC$log2FoldChange
      names(GIC_genelist) <- GIC$ensemble_gene_id
      GIC_genelist <- sort(GIC_genelist,decreasing = TRUE)
      GIC_GSEA <- GSEA(geneList = GIC_genelist,TERM2GENE = GMT, pvalueCutoff = 1, minGSSize = 0, maxGSSize = 500)
      GIC_GSEA$Description
      pathway <- "Tissue_TumorvsNormal_Log2FC>2_p.adjust<0.05"
      gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
      
      pathway <- "Tissue_TumorvsNormal_Log2FC<-2_p.adjust<0.05"
      gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
      
      enrichplot::gseaplot2(GIC_GSEA,c("PBMC_CRCvsHD_Log2FC>0.67_pvalue<0.05","PBMC_CRCvsHD_Log2FC<-0.67_pvalue<0.05",
                                       "Tissue_TumorvsNormal_Log2FC>2.5_padjust<0.05","Tissue_TumorvsNormal_Log2FC<-2.5_padjust<0.05"),
                            base_size = 24,
                            color = c("PBMC_CRCvsHD_Log2FC>0.67_pvalue<0.05"="blue",
                                      "PBMC_CRCvsHD_Log2FC<-0.67_pvalue<0.05"="dark blue",
                                      "Tissue_TumorvsNormal_Log2FC>2.5_padjust<0.05"="red",
                                      "Tissue_TumorvsNormal_Log2FC<-2.5_padjust<0.05"="dark red"),
                            #color = "black",
                            rel_heights = c(2, 0.25, 1),
                            subplots = c(1,2,3),
                            ES_geom = "line",
                            pvalue_table = TRUE)
      library(enrichplot)
      library(gridExtra)
      library(grid)
      trace("gseaplot2", edit = TRUE)
      
      ridgeplot(GIC_GSEA)
    }
  }

#Plot immune geneset on MSigDB
  {
    #prepare GSEA reference geneset (KEGG)
    MsigDB_forGSEA <- read.delim("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/MSigDB/c7.immunesigdb.v7.5.1.symbols.gmt",header  = FALSE,sep = "\t",row.names = 1)
    MsigDB_forGSEA <- MsigDB_forGSEA[,-1]
    MsigDB_forGSEA <- t(MsigDB_forGSEA)
    GMT <- melt(MsigDB_forGSEA)
    GMT <- GMT[,-1]
    colnames(GMT) <- c("DESCRPTION","hgnc_id")
    
    values <- GMT$hgnc_id
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    candidate_list_raw <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), 
                                filters = "hgnc_symbol", 
                                values=unique(values), mart= mart,useCache = FALSE)
    GMT_MsigDB <- left_join(GMT,candidate_list_raw, by = c("hgnc_id"="hgnc_symbol"))
    
    GMT_MsigDB <- GMT_MsigDB[,c("DESCRPTION","ensembl_gene_id")]
    
    #GIC
    {
      GIC_Differential <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/output/Expression_tissue_CRC_TvsN_edger_exact.txt",header = T,sep="\t",quote = "")
      GIC_Differential$ensemble_gene_id <- as.character(lapply(strsplit(as.character(rownames(GIC_Differential)),".",fixed = TRUE),function(x) x[1]))
      GIC_Differential <- aggregate(GIC_Differential[,-which(colnames(GIC_Differential)=="ensemble_gene_id")], list(ensemble_gene_id=GIC_Differential$ensemble_gene_id), FUN = mean)
      GIC <- GIC_Differential[grep("ENSG",GIC_Differential$ensemble_gene_id),c("ensemble_gene_id","log2FoldChange")]
      GIC_genelist <- GIC$log2FoldChange
      names(GIC_genelist) <- GIC$ensemble_gene_id
      GIC_genelist <- sort(GIC_genelist,decreasing = TRUE)
      GIC_GSEA <- GSEA(geneList = GIC_genelist,TERM2GENE = GMT_MsigDB, pvalueCutoff = 1, minGSSize = 20, maxGSSize = 203)
      
      result <- GIC_GSEA@result
      
      #tolower(result[(result$NES > 0) & (result$pvalue < 0.05),]$ID)
      #grep("tcell",tolower(result[(result$NES > 1) & (result$p.adjust < 0.05),]$ID),value = TRUE)
      
      #pathway <- toupper(c("gse46025_wt_vs_foxo1_ko_klrg1_low_cd8_effector_tcell_dn","gse8835_healthy_vs_cll_cd8_tcell_up","gse9650_effector_vs_exhausted_cd8_tcell_up","gse9650_effector_vs_memory_cd8_tcell_dn"))
      #gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
      
      #pathway <- toupper(c("gse46025_wt_vs_foxo1_ko_klrg1_low_cd8_effector_tcell_dn","gse9650_effector_vs_exhausted_cd8_tcell_up","gse9650_effector_vs_memory_cd8_tcell_dn"))
      #gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
      
      #Foxo proteins cooperatively control the differentiation of Foxp3+ regulatory T cells - Foxo proteins regulated the expression of additional Treg cell–associated genes and were essential for inhibiting the acquisition of effector T cell characteristics by Treg cells.
      enrichplot::gseaplot2(GIC_GSEA,toupper(c("gse19825_naive_vs_day3_eff_cd8_tcell_up","gse21678_wt_vs_foxo1_foxo3_ko_treg_up",
                                               "gse9650_effector_vs_exhausted_cd8_tcell_up")),
                            base_size = 24,
                            color = c("GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP"="Blue",
                                      #"GSE8835_HEALTHY_VS_CLL_CD8_TCELL_UP"="dark blue",
                                      "GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP"="red",
                                      "GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP"="dark red"),
                            #color = "black",
                            rel_heights = c(2, 0.25, 1),
                            subplots = c(1,2,3),
                            ES_geom = "line",
                            pvalue_table = FALSE)
      library(enrichplot)
      library(gridExtra)
      library(grid)
      trace("gseaplot2", edit = TRUE)
    }
  
  }
  
#exhaustion immune geneset of 2016-Science-Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq
  {
    #prepare GMT
    GMT <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/gene/T cell Exaustion Consensus genes.csv",header = TRUE)
    GMT_exhaustion <- data.frame("DESCRPTION"="T cell exhaustion","ensembl_gene_id"=GMT$ensembl)
    
    #GIC
    {
      GIC_Differential <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/output/Expression_tissue_CRC_TvsN_edger_exact.txt",header = T,sep="\t",quote = "")
      GIC_Differential$ensemble_gene_id <- as.character(lapply(strsplit(as.character(rownames(GIC_Differential)),".",fixed = TRUE),function(x) x[1]))
      GIC_Differential <- aggregate(GIC_Differential[,-which(colnames(GIC_Differential)=="ensemble_gene_id")], list(ensemble_gene_id=GIC_Differential$ensemble_gene_id), FUN = mean)
      GIC <- GIC_Differential[grep("ENSG",GIC_Differential$ensemble_gene_id),c("ensemble_gene_id","log2FoldChange")]
      GIC_genelist <- GIC$log2FoldChange
      names(GIC_genelist) <- GIC$ensemble_gene_id
      GIC_genelist <- sort(GIC_genelist,decreasing = TRUE)
      GIC_GSEA <- GSEA(geneList = GIC_genelist,TERM2GENE = GMT_exhaustion, pvalueCutoff = 1, minGSSize = 0, maxGSSize = 300)
      
      result <- GIC_GSEA@result
      
      enrichplot::gseaplot2(GIC_GSEA,"T cell exhaustion",
                            base_size = 24,
                            color = c("T cell exhaustion"="Blue"),
                            #color = "black",
                            rel_heights = c(2, 0.25, 1),
                            subplots = c(1,2,3),
                            ES_geom = "line",
                            pvalue_table = FALSE)
      library(enrichplot)
      library(gridExtra)
      library(grid)
      trace("gseaplot2", edit = TRUE)
    }
    
  }
  
#Plot cancer geneset on MSigDB
  {
    #prepare GSEA reference geneset (KEGG)
    MsigDB_forGSEA <- read.delim("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/MSigDB/c6.all.v7.5.1.symbols.gmt",header  = FALSE,sep = "\t",row.names = 1)
    MsigDB_forGSEA <- MsigDB_forGSEA[,-1]
    MsigDB_forGSEA <- t(MsigDB_forGSEA)
    GMT <- melt(MsigDB_forGSEA)
    GMT <- GMT[,-1]
    colnames(GMT) <- c("DESCRPTION","hgnc_id")
    
    values <- GMT$hgnc_id
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    candidate_list_raw <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), 
                                filters = "hgnc_symbol", 
                                values=unique(values), mart= mart,useCache = FALSE)
    GMT_MsigDB <- left_join(GMT,candidate_list_raw, by = c("hgnc_id"="hgnc_symbol"))
    
    GMT_MsigDB <- GMT_MsigDB[,c("DESCRPTION","ensembl_gene_id")]
    
    #GIC
    {
      GIC_Differential <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/output/Expression_GIvsNC_edger_exact.txt",header = T,sep="\t",quote = "")
      GIC_Differential$ensemble_gene_id <- as.character(lapply(strsplit(as.character(rownames(GIC_Differential)),".",fixed = TRUE),function(x) x[1]))
      GIC_Differential <- aggregate(GIC_Differential[,-which(colnames(GIC_Differential)=="ensemble_gene_id")], list(ensemble_gene_id=GIC_Differential$ensemble_gene_id), FUN = mean)
      GIC <- GIC_Differential[grep("ENSG",GIC_Differential$ensemble_gene_id),c("ensemble_gene_id","log2FoldChange")]
      GIC_genelist <- GIC$log2FoldChange
      names(GIC_genelist) <- GIC$ensemble_gene_id
      GIC_genelist <- sort(GIC_genelist,decreasing = TRUE)
      GIC_GSEA <- GSEA(geneList = GIC_genelist,TERM2GENE = GMT_MsigDB, pvalueCutoff = 1, minGSSize = 20, maxGSSize = 203)
      
      result <- GIC_GSEA@result
      
      tolower(result[(result$NES > 0) & (result$pvalue < 0.05),]$ID)
      grep("crc|gastric|stomach|colorectal",tolower(result[(result$NES > 1) & (result$p.adjust < 0.05),]$ID),value = TRUE)
      
      pathway <- toupper(c("gse46025_wt_vs_foxo1_ko_klrg1_low_cd8_effector_tcell_dn","gse8835_healthy_vs_cll_cd8_tcell_up","gse9650_effector_vs_exhausted_cd8_tcell_up","gse9650_effector_vs_memory_cd8_tcell_dn"))
      gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
      
      pathway <- toupper(c("gse46025_wt_vs_foxo1_ko_klrg1_low_cd8_effector_tcell_dn","gse9650_effector_vs_exhausted_cd8_tcell_up","gse9650_effector_vs_memory_cd8_tcell_dn"))
      gseaplot(GIC_GSEA, geneSetID = which(GIC_GSEA$Description==pathway), title = pathway)
      
      #Foxo proteins cooperatively control the differentiation of Foxp3+ regulatory T cells - Foxo proteins regulated the expression of additional Treg cell–associated genes and were essential for inhibiting the acquisition of effector T cell characteristics by Treg cells.
      enrichplot::gseaplot2(GIC_GSEA,toupper(c("gse19825_naive_vs_day3_eff_cd8_tcell_up","gse21678_wt_vs_foxo1_foxo3_ko_treg_up",
                                               "gse9650_effector_vs_exhausted_cd8_tcell_up")),
                            base_size = 24,
                            color = c("GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP"="Blue",
                                      #"GSE8835_HEALTHY_VS_CLL_CD8_TCELL_UP"="dark blue",
                                      "GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP"="red",
                                      "GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP"="dark red"),
                            #color = "black",
                            rel_heights = c(2, 0.25, 1),
                            subplots = c(1,2,3),
                            ES_geom = "line",
                            pvalue_table = FALSE)
      library(enrichplot)
      library(gridExtra)
      library(grid)
      trace("gseaplot2", edit = TRUE)
    }
  }
}

#CPM
{
count_matrix <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/Expression/matrix/20211113_count_matrix.txt",sep = "\t", check.names = FALSE)
count_matrix$`Gene|Length` <- as.character(lapply(strsplit(as.character(count_matrix$`Gene|Length`),".",fixed = TRUE),function(x) x[1]))
count_matrix <- aggregate(count_matrix[,-1], list(Gene=count_matrix[,1]), FUN = sum)
rownames(count_matrix) <- count_matrix$Gene
count_matrix <- count_matrix[,-1]
count_matrix <- count_matrix[,grep("pico",colnames(count_matrix))]
count_matrix <- count_matrix[,-grep("mix..pico",colnames(count_matrix))]
colnames(count_matrix) <- gsub("-",".",fixed = TRUE, colnames(count_matrix))

cpm <- apply(count_matrix, 2, function(x) x / sum(as.numeric(x)) * 10^6)

values <- as.character(lapply(strsplit(rownames(cpm),".",fixed = TRUE),function(x) x[1]))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
annotations <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), 
                     filters = "ensembl_gene_id", 
                     values=unique(values), mart= mart,useCache = FALSE)
annotations$hgnc_symbol <- gsub("-",".",fixed = TRUE, annotations$hgnc_symbol)
write.csv(annotations,"annotations.csv")
annotations <- read.csv("annotations.csv",row.names = 1)
cpm <- as.data.frame(cpm)
cpm$ensembl_gene_id <- rownames(cpm)
cpm_new <- left_join(cpm,annotations,by = c('ensembl_gene_id'='ensembl_gene_id'))

genes <- cpm_new$`hgnc_symbol`
cpm_new <- aggregate(cpm_new[,-c(ncol(cpm_new),ncol(cpm_new)-1)], list(Gene=genes), FUN = sum)
rownames(cpm_new) <- cpm_new$Gene

write.table(cpm_new[,-1],"/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/multiomics_paired_20211113_CPM.txt",sep = "\t", quote = FALSE, row.names = TRUE,col.names = TRUE)
}
#CPM for ipico
{
  count_matrix <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_rawcount/ipico_gencode.cutoff_filter.txt",sep = "\t", check.names = FALSE)
  count_matrix$Gene_name <- as.character(lapply(strsplit(as.character(rownames(count_matrix)),"|",fixed = TRUE),function(x) x[3]))
  count_matrix <- aggregate(count_matrix[,-which(colnames(count_matrix)=="Gene_name")], list(Gene_name=count_matrix$Gene_name), FUN = sum)
  rownames(count_matrix) <- count_matrix$Gene_name
  count_matrix <- count_matrix[,-which(colnames(count_matrix)=="Gene_name")]
  colnames(count_matrix) <- gsub("-","_",fixed = TRUE, colnames(count_matrix))
  
  cpm <- apply(count_matrix, 2, function(x) x / sum(as.numeric(x)) * 10^6)
  
  write.table(cpm,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/ipico_gencode.cutoff_filter_CPM.csv",sep = ",", quote = FALSE, row.names = TRUE,col.names = TRUE)
}

#ActivePathways
{
library(ActivePathways)
fname_scores <- system.file("extdata", "Adenocarcinoma_scores_subset.tsv", 
                            package = "ActivePathways")
fname_GMT = system.file("extdata", "hsapiens_REAC_subset.gmt",
                        package = "ActivePathways")

dat <- as.matrix(read.table(fname_scores, header = TRUE, row.names = 'Gene'))
dat[is.na(dat)] <- 1

ActivePathways(dat, fname_GMT)


#for activepathway
{
Alt <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/15.examples/Alt.promoter/output/AlternativePromoter_GIvsNC_wilcox.txt",header = TRUE, check.names = FALSE)

Alt$Gene <- as.character(lapply(strsplit(as.character(rownames(Alt)),".",fixed = TRUE),function(x) x[1]))

Alt2 <- aggregate(Alt[,-which(colnames(Alt)=="Gene")], list(Gene=Alt$Gene), FUN = min)
rownames(Alt2) <- Alt2$Gene
Alt2 <- Alt2[,-1]

write.table(Alt2,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/15.examples/Alt.promoter/output/AlternativePromoter_GIvsNC_wilcox_min.txt",row.names = TRUE,quote = FALSE)
}

}

#TIDE-calculate immune disfunctional score: http://tide.dfci.harvard.edu/
#prepare matrix with entrezgene_id and split by tab
#log2(TPM+1) value and used all sample as normalization control
{
  #raw_matrix <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Tissue_TPM.txt",row.names = 1,sep = "\t",header = TRUE,check.names = FALSE)
  raw_matrix <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/GSE174302_intron-spanning-available-samples_TPM.txt",row.names = 1,sep = "\t",header = TRUE,check.names = FALSE)
  raw_matrix <- raw_matrix[grep("ENSG",rownames(raw_matrix)),]
  
  #raw_matrix$ensembl_id <- as.character(lapply(strsplit(as.character(rownames(raw_matrix)),".",fixed = TRUE),function(x) x[1])) 
  raw_matrix$ensembl_id <- as.character(lapply(strsplit(as.character(rownames(raw_matrix)),"|",fixed = TRUE),function(x) x[1])) 
  
  library(biomaRt)
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", mirror = "asia")
  gene_names <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                      filters = "ensembl_gene_id",
                      values=raw_matrix$ensembl_id, mart= mart,useCache = FALSE)
  raw_matrix <- left_join(raw_matrix,gene_names,by = c("ensembl_id"="ensembl_gene_id"))
  
  raw_matrix$entrezgene_id <- as.character(raw_matrix$entrezgene_id)
  raw_matrix[raw_matrix==""] <- NA
  raw_matrix <- raw_matrix %>%
    mutate(entrezgene_id = coalesce(entrezgene_id,ensembl_id))
  
  raw_matrix <- raw_matrix[-grep("ENSG",raw_matrix$entrezgene_id),]
  raw_matrix <- raw_matrix[,-which(colnames(raw_matrix)=="ensembl_id")]

  aggregated_matrix <- aggregate(raw_matrix[,-which(colnames(raw_matrix)=="entrezgene_id")], list(Gene=raw_matrix[,which(colnames(raw_matrix)=="entrezgene_id")]), FUN = sum)
  rownames(aggregated_matrix) <- aggregated_matrix$Gene
  aggregated_matrix <- aggregated_matrix[,-which(colnames(aggregated_matrix)=="Gene")]
  aggregated_matrix <- log2(aggregated_matrix+1)
  
  test <- apply(aggregated_matrix,1,function(x) x/mean(x))
  
  write.table(as.data.frame(t(test)),"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/GSE174302_matrix_forTIDE.txt",sep = "\t",quote = FALSE)
  
  #then: tidepy -o GSE174302_TIDE.txt -c Other GSE174302_matrix_forTIDE.txt
}


##stromal cell fraction by epic (elife)
#devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
{
library("EPIC")
#raw_matrix <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Plasma_TPM.txt",row.names = 1,sep = "\t",header = TRUE,check.names = FALSE)
raw_matrix <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/GSE174302_intron-spanning-available-samples_TPM.txt",row.names = 1,sep = "\t",header = TRUE,check.names = FALSE)
raw_matrix <- raw_matrix[grep("ENSG",rownames(raw_matrix)),]

#raw_matrix$ensembl_id <- as.character(lapply(strsplit(as.character(rownames(raw_matrix)),".",fixed = TRUE),function(x) x[1])) 
raw_matrix$ensembl_id <- as.character(lapply(strsplit(as.character(rownames(raw_matrix)),"|",fixed = TRUE),function(x) x[1])) 

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", mirror = "asia")
gene_names <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values=raw_matrix$ensembl_id, mart= mart,useCache = FALSE)
raw_matrix <- left_join(raw_matrix,gene_names,by = c("ensembl_id"="ensembl_gene_id"))

raw_matrix$hgnc_symbol <- as.character(raw_matrix$hgnc_symbol)
raw_matrix[raw_matrix==""] <- NA
raw_matrix <- raw_matrix %>%
  mutate(hgnc_symbol = coalesce(hgnc_symbol,ensembl_id))

raw_matrix <- raw_matrix[-grep("ENSG",raw_matrix$hgnc_symbol),]
raw_matrix <- raw_matrix[,-which(colnames(raw_matrix)=="ensembl_id")]

aggregated_matrix <- aggregate(raw_matrix[,-which(colnames(raw_matrix)=="hgnc_symbol")], list(Gene=raw_matrix[,which(colnames(raw_matrix)=="hgnc_symbol")]), FUN = sum)
rownames(aggregated_matrix) <- aggregated_matrix$Gene
aggregated_matrix <- aggregated_matrix[,-which(colnames(aggregated_matrix)=="Gene")]

result <- EPIC(bulk = aggregated_matrix, reference = "TRef")

write.table(result$mRNAProportions,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/GSE174302_EPIC_mRNAproportion_TRef.txt",sep = "\t",quote = FALSE)
write.table(result$cellFractions,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/GSE174302_EPIC_cellFraction_TRef.txt",sep = "\t",quote = FALSE)
write.table(result$fit.gof,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/GSE174302_EPIC_fitquality_TRef.txt",sep = "\t",quote = FALSE)
}

#immune cell fraction barplot
#plasma
{
  #read in immune fractions
  {
  Immune_fraction <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/res_cibersort_multiomics_plasma_PBMC_tissue_LM22_abs_20220608_perm1000.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
  Immune_fraction2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/res_cibersort_GSE174302_LM22_abs_20220614_perm1000.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
  Immune_fraction <- rbind(Immune_fraction,Immune_fraction2)
  
  Immune_fraction <- Immune_fraction[,-grep("P-value|Correlation|RMSE|Absolute score",colnames(Immune_fraction))]
  {
    Immune_fraction$`B cells` <- Immune_fraction$`B cells naive`+Immune_fraction$`B cells memory`
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="B cells naive")]
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="B cells memory")]
    Immune_fraction$`T cells CD4` <- Immune_fraction$`T cells CD4 naive`+Immune_fraction$`T cells CD4 memory resting`+Immune_fraction$`T cells CD4 memory activated`+Immune_fraction$`T cells follicular helper`+Immune_fraction$`T cells regulatory (Tregs)`
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="T cells CD4 naive")]
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="T cells CD4 memory resting")]
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="T cells CD4 memory activated")]
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="T cells follicular helper")]
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="T cells regulatory (Tregs)")]
    Immune_fraction$`NK cells` <- Immune_fraction$`NK cells resting`+Immune_fraction$`NK cells activated`
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="NK cells resting")]
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="NK cells activated")]
    Immune_fraction$Macrophages <- Immune_fraction$`Macrophages M0`+Immune_fraction$`Macrophages M1`+Immune_fraction$`Macrophages M2`
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Macrophages M0")]
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Macrophages M1")]
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Macrophages M2")]
    Immune_fraction$`Dendritic cells` <- Immune_fraction$`Dendritic cells resting`+Immune_fraction$`Dendritic cells activated`
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Dendritic cells resting")]
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Dendritic cells activated")]
    Immune_fraction$`Mast cells` <- Immune_fraction$`Mast cells resting`+Immune_fraction$`Mast cells activated`
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Mast cells resting")]
    Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Mast cells activated")]
  }
  #Immune_fraction <- Immune_fraction[-grep("CRC.PKU.15|STAD.PKU.23|STAD.PKU.36",rownames(Immune_fraction)),]
  #Immune_fraction <- Immune_fraction[grep("STAD",rownames(Immune_fraction)),]
  Immune_fraction <- Immune_fraction[-which(rowSums(Immune_fraction)==0),] # This step will remove samples that have no immune cells fraction
  Immune_fraction <- Immune_fraction/rowSums(Immune_fraction)
  Immune_fraction$ID <- rownames(Immune_fraction)
  Immune_fraction$ID <- gsub(".","-",fixed = TRUE,Immune_fraction$ID)
  rownames(Immune_fraction) <- gsub(".","-",fixed = TRUE,rownames(Immune_fraction))
  }
  
  rank1 <- Immune_fraction[Immune_fraction$`T cells gamma delta`>0.0,]
  rank2 <- Immune_fraction[Immune_fraction$`T cells gamma delta`<=0.0,]
  final_rank <- c(rank1[order(rank1$`T cells gamma delta`+rank1$`B cells`+rank1$`Plasma cells`,decreasing = TRUE),]$ID,
                  rank2[order(rank2$`B cells`+rank2$`Plasma cells`,decreasing = TRUE),]$ID)
  final_rank <- Immune_fraction[order(Immune_fraction$`T cells gamma delta`,Immune_fraction$`T cells CD4`,decreasing = TRUE),]$ID
  
  Immune_fraction_forplot <- melt(Immune_fraction,id.vars = "ID")
  Immune_fraction_forplot <- Immune_fraction_forplot[-grep("NC|PBMC|-T|-N|HCC|LUAD|ESCA|Tumor|Normal|mix|CRC-PKU-5|NC-PKU-mix17|STAD-PKU-4",Immune_fraction_forplot$ID),]
  #Immune_fraction_forplot <- Immune_fraction_forplot[-grep("CRC-PKU-15|STAD-PKU-23|STAD-PKU-36|STAD-PKU-24|STAD-PKU-35",Immune_fraction_forplot$ID),]
  
  Immune_fraction_forplot$ID <- factor(Immune_fraction_forplot$ID,levels = final_rank)
  #Immune_fraction_forplot$ID <- factor(Immune_fraction_forplot$ID,levels = Immune_fraction[order((Immune_fraction$`T cells gamma delta`-Immune_fraction$`CD8 T cells`),decreasing = FALSE),]$ID)
  
  
  #read in clinical data
  {
    clinical <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Sample_clinical_information_x.csv",header = TRUE,check.names = FALSE,row.names = 1)
    clinical$ID <- paste0(rownames(clinical),"-pico")
    clinical <- clinical[-grep("NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",clinical$ID),]
    
    clinical2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Sample_clinical_information_GSE174302.csv",header = TRUE,check.names = FALSE,row.names = 1)
    rownames(clinical2) <- gsub("CRC-2410911","CRC-2410966",rownames(clinical2))
    clinical2$ID <- rownames(clinical2)
    clinical <- rbind(clinical,clinical2)
    
    
    #clinical <- clinical[-grep("NC|PBMC|Tumor|Normal|mix",clinical$ID),]
    rownames(clinical) <- clinical$ID
    clinical <- clinical[as.character(unique(Immune_fraction_forplot$ID)),]
    clinical$ID <- factor(clinical$ID,levels=levels(Immune_fraction_forplot$ID))
    #clinical$ID <- factor(clinical$ID,levels=c(clinical[order(clinical$`Tumor size`,decreasing = TRUE),]$ID))
    #clinical$ID <- factor(clinical$ID,levels=c(clinical[order(clinical$T,decreasing = TRUE),]$ID))
    #Immune_fraction_forplot$ID <- factor(Immune_fraction_forplot$ID,levels = levels(clinical$ID))
  }
  
  #read in TIDE scores
  { 
    TIDE <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Plasma_TIDE.txt",sep="\t",header = TRUE,check.names = FALSE,row.names = 1)
    TIDE$ID <- rownames(TIDE)
    TIDE <- TIDE[-grep("NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",TIDE$ID),]
    
    TIDE2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/GSE174302_TIDE.txt",sep="\t",header = TRUE,check.names = FALSE,row.names = 1)
    TIDE2$ID <- rownames(TIDE2)
    
    TIDE <- rbind(TIDE,TIDE2)
    
    TIDE$group <- as.character(lapply(strsplit(TIDE$ID,"-"), function(x) x[1]))
    TIDE$group <- gsub("NC","HD",TIDE$group)
    
    
    #my_comparisons <- list(c("STAD","HD"),c("CRC","HD"))
    #Tumor_comparison <- 
    #ggplot(TIDE,aes(x=group,y=Dysfunction,fill = group))+
    #  geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
    #  geom_point(size = 1, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
    #  scale_fill_manual(values=c("STAD"="red","CRC"="#FCB514","HD"="blue")) +
    #  theme_bw()+
    #  ylab("Dysfunction score")+
    #  xlab("")+
    #  theme(#legend.position="right",
    #    panel.grid=element_blank(),
    #    panel.border=element_blank(),
    #    axis.line = element_line(size=1, colour = "black"),
    #    legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
    #    legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
    #    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    #    axis.text.x = element_blank(),
    #    #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
    #    axis.text.y = element_text(face="bold",  color="black", size=24),
    #    axis.title.x = element_text(face="bold", color="black", size=24),
    #    axis.title.y = element_text(face="bold",color="black", size=24))+
    #  stat_compare_means(comparisons = my_comparisons,
    #                     method = "wilcox.test",
    #                     size = 13,
    #                     vjust = 0.6,
    #                     method.args = list(alternative = "two.sided",paired = TRUE),
    #                     label = "p.signif")
    #ggsave(plot = Tumor_comparison, filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Tumor_dysfunction_comparison.pdf",device = "pdf",width = 3.8,height = 5)
    
    TIDE <- TIDE[as.character(unique(Immune_fraction_forplot$ID)),]
    TIDE$group <- factor(TIDE$group,levels = c("CRC","STAD","HD"))
    TIDE$ID <- factor(TIDE$ID,levels=levels(Immune_fraction_forplot$ID))
  }
  
  #read in EPIC scores
  {
    EPIC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Plasma_EPIC_cellFraction_TRef.txt",sep="\t",header = TRUE,check.names = FALSE,row.names = 1)
    EPIC$ID <- rownames(EPIC)
    EPIC <- EPIC[-grep("NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",EPIC$ID),]
    
    EPIC2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/GSE174302_EPIC_cellFraction_TRef.txt",sep="\t",header = TRUE,check.names = FALSE,row.names = 1)
    EPIC2$ID <- rownames(EPIC2)
    EPIC <- rbind(EPIC,EPIC2)
    
    EPIC <- EPIC[as.character(unique(Immune_fraction_forplot$ID)),]
    EPIC$ID <- factor(EPIC$ID,levels = levels(Immune_fraction_forplot$ID))
    colnames(EPIC) <- paste0(colnames(EPIC),"_EPIC")
  }
  
  #read in gene expression signatures
  {
    T_cell_receptor_signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/T cell receptor signaling pathway_matrix/Multiomics_20211113/Multiomics_20211113_T cell receptor signaling pathway_log2sum.txt",sep = "\t")
    T_cell_receptor_signature2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/T cell receptor signaling pathway_matrix/GSE174302/GSE174302_T cell receptor signaling pathway_log2sum.txt",sep = "\t")
    T_cell_receptor_signature <- rbind(T_cell_receptor_signature,T_cell_receptor_signature2)
    
    B_cell_receptor_signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/B cell receptor signaling pathway_matrix/Multiomics_20211113/Multiomics_20211113_B cell receptor signaling pathway_log2sum.txt",sep = "\t")
    B_cell_receptor_signature2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/B cell receptor signaling pathway_matrix/GSE174302/GSE174302_B cell receptor signaling pathway_log2sum.txt",sep = "\t")
    B_cell_receptor_signature <- rbind(B_cell_receptor_signature,B_cell_receptor_signature2)
    
    Naive_T_cell_signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP_matrix/Multiomics_20211113/Multiomics_20211113_GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP_log2sum.txt",sep = "\t")
    Naive_T_cell_signature2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP_matrix/GSE174302/GSE174302_GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP_log2sum.txt",sep = "\t")
    Naive_T_cell_signature <- rbind(Naive_T_cell_signature, Naive_T_cell_signature2)
    
    FOXO_regulated_signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP_matrix/Multiomics_20211113/Multiomics_20211113_GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP_log2sum.txt",sep = "\t")
    FOXO_regulated_signature2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP_matrix/GSE174302/GSE174302_GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP_log2sum.txt",sep = "\t")
    FOXO_regulated_signature <- rbind(FOXO_regulated_signature,FOXO_regulated_signature2)
    
    Effector_T_cell_receptor_signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP_matrix/Multiomics_20211113/Multiomics_20211113_GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP_log2sum.txt",sep = "\t")
    Effector_T_cell_receptor_signature2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP_matrix/GSE174302/GSE174302_GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP_log2sum.txt",sep = "\t")
    Effector_T_cell_receptor_signature <- rbind(Effector_T_cell_receptor_signature,Effector_T_cell_receptor_signature2)
    
    Gene_signature <- data.frame(row.names = as.character(unique(Immune_fraction_forplot$ID)),
    "T_cell_receptor_signature"=T_cell_receptor_signature[as.character(unique(Immune_fraction_forplot$ID)),],
    "B_cell_receptor_signature"=B_cell_receptor_signature[as.character(unique(Immune_fraction_forplot$ID)),],
    "Naive_T_cell_signature"=Naive_T_cell_signature[unique(as.character(Immune_fraction_forplot$ID)),],
    "FOXO_regulated_signature"=FOXO_regulated_signature[unique(as.character(Immune_fraction_forplot$ID)),],
    "Effector_T_cell_receptor_signature"=Effector_T_cell_receptor_signature[unique(as.character(Immune_fraction_forplot$ID)),]
    )
  }
  
  #Correlation and comparison
  {
  forCorrlation <- cbind(Immune_fraction[unique(as.character(Immune_fraction_forplot$ID)),],TIDE,EPIC,Gene_signature,clinical)

  forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="Stage")]
  forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="ID")]
  forCorrlation <- forCorrlation[,-grep("ID_EPIC|ID.1|ID.2|No benefits|MSI Score|Patient ID|MSI/MMR|group|Diagnosis|CTL.flag|Responder",colnames(forCorrlation))]
  #CRC
  forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="M")]
  forCorrlation <- forCorrlation[,-grep("No benefits|MSI Score|Type|Patient ID|MSI/MMR|group|Diagnosis|CTL.flag|Responder|Multiple primary|MUC2|MUC5AC|EBER|EGFR|CDNA2",colnames(forCorrlation))]
  #STAD
  #forCorrlation <- forCorrlation[,-grep("No benefits|MSI Score|Type|Patient ID|MSI/MMR|group|Diagnosis|CTL.flag|Responder|Right half|EGFR|MLH1|MSH2|MSH6|PMS2",colnames(forCorrlation))]
  
  
  {
  forCorrlation$Stage_simplified <- factor(forCorrlation$Stage_simplified,levels = c("early","late","x"))
  #forCorrlation$`Multiple primary cancer(Yes/No)` <- factor(forCorrlation$`Multiple primary cancer(Yes/No)`,levels = c())
  #forCorrlation$`Right half / left half` <- factor(forCorrlation$`Right half / left half`,levels = c())
  #forCorrlation$`Tumor size` <- factor(forCorrlation$,levels = c())
  forCorrlation$`T` <- factor(forCorrlation$`T`,levels = c("Tis","1","1a","1b","2","3","4","4a","4b","x"))
  #forCorrlation$`N` <- factor(forCorrlation$,levels = c())
  #forCorrlation$`M` <- factor(forCorrlation$,levels = c())
  forCorrlation$`Vascular tumor thrombus` <- factor(forCorrlation$`Vascular tumor thrombus`,levels = c("No","Yes","x"))
  forCorrlation$`Neurological invasion` <- factor(forCorrlation$`Neurological invasion`,levels = c("No","Yes","x"))
  forCorrlation$`Tumor deposition/cancer nodules` <- factor(forCorrlation$`Tumor deposition/cancer nodules`,levels = c("No","Yes","x"))
  #forCorrlation$`CEA(HE)` <- factor(forCorrlation$,levels = c())
  #forCorrlation$`CA199(HE)` <- factor(forCorrlation$,levels = c())
  #forCorrlation$`EGFR` <- factor(forCorrlation$,levels = c())
  #forCorrlation$`HER2` <- factor(forCorrlation$,levels = c())
  #forCorrlation$`CDx2` <- factor(forCorrlation$,levels = c())
  #forCorrlation$P53 <- factor(forCorrlation$,levels = c())
  #forCorrlation$`TOP IIa` <- factor(forCorrlation$,levels = c())
  #forCorrlation$Ki67 <- factor(forCorrlation$,levels = c())
  #forCorrlation$`MMR(p/d)` <- factor(forCorrlation$,levels = c())
  #forCorrlation$MLH1 <- factor(forCorrlation$,levels = c())
  #forCorrlation$PMS2 <- factor(forCorrlation$,levels = c())
  #forCorrlation$MSH2 <- factor(forCorrlation$,levels = c())
  #forCorrlation$MSH6 <- factor(forCorrlation$,levels = c())
  #forCorrlation$MUC2 <- factor(forCorrlation$,levels = c())
  #forCorrlation$MUC5AC <- factor(forCorrlation$,levels = c())
  #forCorrlation$EBER <- factor(forCorrlation$,levels = c())
  }
  
  #method1
  {
  correlation_method <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation/GSE174302_STAD_"
  #enmuerate_correlation(forCorrlation,correlation_method,output_dir)
  forCorrlation <- forCorrlation[-grep("pico",rownames(forCorrlation)),]
  enmuerate_correlation(forCorrlation[grep("STAD",rownames(forCorrlation)),],correlation_method,output_dir)
  
  result <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation/CRC_STAD_Result_final.csv",header = TRUE, row.names = 1,check.names = FALSE)
  result$column <- as.character(lapply(strsplit(as.character(result$DataNames)," vs ",fixed = TRUE), function(x) x[1]))
  result$row <- as.character(lapply(strsplit(as.character(result$DataNames)," vs ",fixed = TRUE), function(x) x[2]))
  R <- data.frame(matrix(nrow=length(colnames(forCorrlation)),ncol = length(colnames(forCorrlation))))
  rownames(R) = colnames(forCorrlation)
  colnames(R) = colnames(forCorrlation)
  pvalue <- data.frame(matrix(nrow=length(colnames(forCorrlation)),ncol = length(colnames(forCorrlation))))
  rownames(pvalue) = colnames(forCorrlation)
  colnames(pvalue) = colnames(forCorrlation)
  
  i=1
  while(i<=nrow(result)){
    pvalue[result[i,"column"],result[i,"row"]] <- as.numeric(as.character(result[i,"Pvalue"]))
    pvalue[result[i,"row"],result[i,"column"]] <- as.numeric(as.character(result[i,"Pvalue"]))
    R[result[i,"column"],result[i,"row"]] <- as.numeric(as.character(result[i,"R"]))
    R[result[i,"row"],result[i,"column"]] <- as.numeric(as.character(result[i,"R"]))
    i=i+1
  }
  
  i=1
  while(i<=nrow(R)){
    pvalue[i,i] <- 0
    R[i,i] <- 1
    i=i+1
  }
  
  R[is.na(R)] <- 0
  pvalue[is.na(pvalue)] <- 1
  
  r = rbind(c('Plasma cells','Type','Mast cells','EBER'),
            c('TIDE','Type','TAM M2','EBER'),
            c('Bcells_EPIC','Type','otherCells_EPIC','EBER'),
            c('T_cell_receptor_signature','Type','Effector_T_cell_receptor_signature','EBER'),
            c('Plasma cells','Type','Effector_T_cell_receptor_signature','Age'),
            c('Plasma cells','Stage_raw','Effector_T_cell_receptor_signature','Multiple primary cancer'),
            c('Plasma cells','Hgb','Effector_T_cell_receptor_signature','PT'),
            c('Plasma cells','CEA','Effector_T_cell_receptor_signature','CA199(HE)'),
            c('Plasma cells','EGFR','Effector_T_cell_receptor_signature','EBER')
  )
  corrplot(corr = as.matrix(R),tl.col="black",order="original",tl.pos = "ld",tl.cex=0.7,tl.srt = 45, type="lower", col = col2(200), #mar = c(1, 1, 1, 1),
           p.mat = as.matrix(pvalue), sig.level = 0.05,insig = "blank",pch.cex = 2) %>% corrRect(namesMat = r)
  }
  
  #method2
  {
  library(corrplot)
  library(RColorBrewer)
  library("Hmisc")
  library(magrittr)
  col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                 "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                 "#4393C3", "#2166AC", "#053061")))
  res2 <- rcorr(as.matrix(sapply(forCorrlation, as.numeric))  ,type="spearman")
  
  res2$r[is.na(res2$r)] <- 1
  res2$P[is.na(res2$P)] <- 0
  
  r = rbind(c('Plasma cells','Type','Mast cells','EBER'),
            c('TIDE','Type','TAM M2','EBER'),
            c('Bcells','Type','otherCells','EBER'),
            c('T_cell_receptor_signature','Type','Effector_T_cell_receptor_signature','EBER'),
            c('Plasma cells','Type','Effector_T_cell_receptor_signature','Age'),
            c('Plasma cells','Stage_raw','Effector_T_cell_receptor_signature','Multiple primary cancer'),
            c('Plasma cells','Hgb','Effector_T_cell_receptor_signature','PT'),
            c('Plasma cells','CEA','Effector_T_cell_receptor_signature','CA199(HE)'),
            c('Plasma cells','EGFR','Effector_T_cell_receptor_signature','EBER')
            )
  
  r = rbind(c('B cells','Gender','Neutrophils','EBER'),
            c('TIDE','Gender','TAM M2','EBER'),
            c('Bcells','Gender','otherCells','EBER'),
            c('T_cell_receptor_signature','Gender','Effector_T_cell_receptor_signature','EBER'),
            c('B cells','Gender','Effector_T_cell_receptor_signature','Age'),
            c('B cells','Stage_raw','Effector_T_cell_receptor_signature','Multiple primary cancer'),
            c('B cells','Hgb','Effector_T_cell_receptor_signature','PT'),
            c('B cells','CEA','Effector_T_cell_receptor_signature','CA199(HE)'),
            c('B cells','HER2','Effector_T_cell_receptor_signature','EBER')
  )
  
  r = rbind(c('B cells','Gender','Neutrophils','MSH6'),
            c('TIDE','Gender','TAM M2','MSH6'),
            c('Bcells','Gender','otherCells','MSH6'),
            c('T_cell_receptor_signature','Gender','Effector_T_cell_receptor_signature','MSH6'),
            c('B cells','Gender','Effector_T_cell_receptor_signature','Age'),
            c('B cells','Stage_raw','Effector_T_cell_receptor_signature','Tumor deposition/cancer nodules'),
            c('B cells','Hgb','Effector_T_cell_receptor_signature','PT'),
            c('B cells','CEA','Effector_T_cell_receptor_signature','CA199(HE)'),
            c('B cells','EGFR','Effector_T_cell_receptor_signature','MSH6')
  )
  
  corrplot(corr = res2$r,tl.col="black",order="original",tl.pos = "ld",tl.cex=0.7,tl.srt = 45, type="lower", col = col2(200), #mar = c(1, 1, 1, 1),
           p.mat = res2$P, sig.level = 0.05,insig = "blank",pch.cex = 2) %>% corrRect(namesMat = r)
  }
  
  #comparison
  {
  forCorrlation$boxplot_group <- NA
  forCorrlation[forCorrlation$`T cells gamma delta`>0,]$boxplot_group <- "gamma_delta_postive"
  forCorrlation[forCorrlation$`T cells gamma delta`<=0,]$boxplot_group <- "gamma_delta_negative"
  
  forCorrlation$Stage_group <- NA
  forCorrlation[forCorrlation$Stage_raw=="1",]$Stage_group <- "Stage I"
  forCorrlation[forCorrlation$Stage_raw=="1A",]$Stage_group <- "Stage I"
  forCorrlation[forCorrlation$Stage_raw=="1B",]$Stage_group <- "Stage I"
  forCorrlation[forCorrlation$Stage_raw=="2A",]$Stage_group <- "Stage II"
  forCorrlation[forCorrlation$Stage_raw=="2B",]$Stage_group <- "Stage II"
  forCorrlation[forCorrlation$Stage_raw=="2C",]$Stage_group <- "Stage II"
  forCorrlation[forCorrlation$Stage_raw=="3A",]$Stage_group <- "Stage III/IV"
  forCorrlation[forCorrlation$Stage_raw=="3B",]$Stage_group <- "Stage III/IV"
  forCorrlation[forCorrlation$Stage_raw=="3C",]$Stage_group <- "Stage III/IV"
  forCorrlation[forCorrlation$Stage_raw=="4",]$Stage_group <- "Stage III/IV"
  
  
  my_comparisons <- list(c("gamma_delta_postive","gamma_delta_negative"))
  my_comparisons <- list(c("early","late"))
  #my_comparisons <- list(c("Stage I","Stage II"),c("Stage I","Stage III"),c("Stage I","Stage IV"))
  my_comparisons <- list(c("Stage I","Stage II"),c("Stage I","Stage III/IV"))
  forCorrlation_boxplot <- forCorrlation[grep("Stomach Cancer",forCorrlation$Type),]
  Tumor_comparison <- 
    ggplot(forCorrlation_boxplot[forCorrlation_boxplot$Stage_simplified!="x",],aes(x=Stage_group, y = CAFs_EPIC, fill = Stage_group))+
    geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
    geom_point(size = 1, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
    scale_fill_manual(values=c("gamma_delta_postive"="dark green","gamma_delta_negative"="white",
                               "early"="sky blue","late"="black",
                               "Stage I"="#90E0EF","Stage II"="#00B4D8","Stage III"="#03045E","Stage III/IV"="#03045E","Stage IV"="black")) +
    theme_bw()+
    ylab("CAFs")+
    xlab("")+
    theme(#legend.position="right",
      #element_text(family = "Arial"),
      legend.position="none",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(size=0.5, colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_blank(),
      #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
      axis.text.y = element_text(color="black", size=24,family = "Arial"),
      axis.title.x = element_text(color="black", size=24,family = "Arial"),
      axis.title.y = element_text(color="black", size=24,family = "Arial"))+
    stat_compare_means(comparisons = my_comparisons,
                       method = "wilcox.test",
                       size = 13,
                       vjust = 0.6,
                       method.args = list(alternative = "less",paired = FALSE),
                       label = "p.signif")
  Tumor_comparison
  ggsave(plot = Tumor_comparison, filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Stage_CAFs_STAD_multiomics+GSE174302_comparison.pdf",device = "pdf",width = 4.5,height = 5)
  }
  }
  
  #plot heatmap and barplot
  {
  Immune_fraction_forplot_selected <- Immune_fraction_forplot[grep("CRC",Immune_fraction_forplot$ID),]
  clinical_selected <- clinical[grep("CRC",clinical$ID),]
  TIDE_selected <- TIDE[grep("CRC",TIDE$ID),]
  EPIC_selected <- EPIC[grep("CRC",EPIC$ID_EPIC),]
  
   Immune_fraction_forplot$variable <- factor(Immune_fraction_forplot$variable, levels = rev(c("T cells gamma delta","T cells CD4","B cells","Plasma cells","T cells CD8",
                                                                                          "NK cells","Monocytes","Macrophages","Dendritic cells",
                                                                                      "Mast cells","Eosinophils","Neutrophils")))

    simpsons <- c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF",
                  "#197EC0FF","#F05C3BFF","#46732EFF","#71D0F5FF","#370335FF","#075149FF",
                  "#C80813FF","#91331FFF","#1A9993FF","#FD8CC1FF")
    
    simpsons <- c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF",
                  "#197EC0FF",alpha("purple",alpha=0.4),alpha("brown",alpha=0.8),"#71D0F5FF","#370335FF","dark green",
                  "#C80813FF","#91331FFF","#1A9993FF","#FD8CC1FF")
  bar <- 
    ggplot(Immune_fraction_forplot_selected,aes(x=ID,y=value,fill=variable))+
    geom_bar(stat = "identity",position = "stack")+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,10,-20,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      #axis.text.x = element_text( color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
      axis.text.y = element_text(color="black", size=16, angle = 90,hjust=0.5, family = "Arial"),
      axis.title.x = element_text(face="bold", color="black", size=20, family = "Arial"),
      axis.title.y = element_text(face="bold",color="black", size=20, family = "Arial"))+
    scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0),limits = c(0,1))+
      scale_fill_manual(values = simpsons)
    #scale_fill_rickandmorty(alpha = 1)
  
  
  
    clinical_TumorSize <-
    ggplot(clinical_selected,aes(x=ID,y=0,fill= as.numeric(as.character(clinical_selected$`Tumor size`))))+
    geom_tile()+
    scale_fill_gradient2(low = "white",
                        mid = "white",
                        high = "red",
                        midpoint = 6)+
    #scale_fill_manual(values=c("-"="white","early"="sky blue","late"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_PLT <-
    ggplot(clinical_selected,aes(x=ID,y=0,fill= as.numeric(as.character(clinical_selected$PLT))))+#`Tumor size`))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "orange")+
    #scale_fill_manual(values=c("-"="white","early"="sky blue","late"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_P53 <-
    ggplot(clinical_selected,aes(x=ID,y=0,fill= P53))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","-"="white","+"="white","++"="dark grey","+++"="black"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_MMR <-
    ggplot(clinical_selected,aes(x=ID,y=0,fill= `MSI/MMR`))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","x"="grey","Normal"="white","Unstable; lack"="#F05C3BFF"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_HER2 <-
    ggplot(clinical_selected,aes(x=ID,y=0,fill= HER2))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","-"="white","+"="white","++"="red","+++"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_Right <-
    ggplot(clinical_selected,aes(x=ID,y=0,fill= `Right half / left half`))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","L"="white","R"="dark red","L+R"="dark red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_Stage <-
    ggplot(clinical_selected,aes(x=ID,y=0,fill= Stage))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","x"="grey","1"="#90E0EF","1A"="#90E0EF","1B"="#90E0EF","2A"="#00B4D8","2B"="#00B4D8","2C"="#00B4D8","3A"="#03045E","3B"="#03045E","3C"="#03045E","4"="black"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_T <-
    ggplot(clinical_selected,aes(x=ID,y=0,fill= `T`))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","Tis"="#CAF0F8","1"="#90E0EF","1a"="#90E0EF","1b"="#00B4D8","2"="#00B4D8","3"="#03045E","4a"="black","4b"="black"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_N <-
    ggplot(clinical_selected,aes(x=ID,y=0,fill= N))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("x"="grey","0"="#CAF0F8","1"="#90E0EF","1a"="#90E0EF","1b"="#90E0EF","2"="#00B4D8","2a"="#00B4D8","2b"="#00B4D8","3a"="#03045E","3b"="#03045E"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_M <-
    ggplot(clinical_selected,aes(x=ID,y=0,fill= as.factor(clinical_selected$M)))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("0"="white","1"="dark red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())

  clinical_Cancer <-
    ggplot(clinical_selected,aes(x=ID,y=0,fill= Type))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("Colorectal Cancer"="#FCB514","Stomach Cancer"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())

  TIDE_Responder <-
    ggplot(TIDE_selected,aes(x=ID,y=0,fill= Responder))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "blue")+
    scale_fill_manual(values=c("True"="dark blue","False"="white"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  TIDE_IFNG <-
    ggplot(TIDE_selected,aes(x=ID,y=0,fill= IFNG))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "sky blue")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  TIDE_Dysfunction <-
    ggplot(TIDE_selected,aes(x=ID,y=0,fill= Dysfunction))+
    geom_tile()+
    scale_fill_gradient2(low="white",mid="white",midpoint = -0.05,high = "dark green")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    #title("MDSC+CAF")+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  TIDE_Exclusion <-
    ggplot(TIDE_selected,aes(x=ID,y=0,fill= Exclusion))+
    geom_tile()+
    scale_fill_gradient2(low="white",mid="white",midpoint = 0,high = "orange")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    #title("MDSC+CAF")+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  TIDE_CAF <-
    ggplot(TIDE_selected,aes(x=ID,y=0,fill= CAF))+
    geom_tile()+
    scale_fill_gradient2(low="white",mid="white",midpoint = 0,high = "dark orange")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    #title("MDSC+CAF")+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  TIDE_TIDE <-
    ggplot(TIDE_selected,aes(x=ID,y=0,fill= TIDE))+
    geom_tile()+
    scale_fill_gradient2(low="white",mid="white",midpoint = 0,high = "orange")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    #title("MDSC+CAF")+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  EPIC_CAFsEndothelial <-
    ggplot(EPIC_selected,aes(x=ID_EPIC,y=0,fill= EPIC_selected$CAFs_EPIC+EPIC_selected$Endothelial_EPIC))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "blue")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  EPIC_CAFs <-
    ggplot(EPIC_selected,aes(x=ID_EPIC,y=0,fill= CAFs_EPIC))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "blue")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  EPIC_Macrophages <-
    ggplot(EPIC_selected,aes(x=ID_EPIC,y=0,fill= Macrophages_EPIC))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "blue")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  EPIC_CAFs <-
    ggplot(EPIC_selected,aes(x=ID_EPIC,y=0,fill= CAFs_EPIC))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "blue")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  EPIC_Bcells <-
    ggplot(EPIC_selected,aes(x=ID_EPIC,y=0,fill= Bcells_EPIC))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "blue")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,0,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  ggarrange(bar,clinical_TumorSize,clinical_PLT,clinical_P53,clinical_MMR,clinical_HER2,clinical_Right,clinical_T,clinical_N,clinical_M,TIDE_Responder,TIDE_IFNG,TIDE_Dysfunction,TIDE_Exclusion,EPIC_CAFsEndothelial,
            ncol = 1,nrow = 14,align = "v",heights = c(100,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5),legend = "none")
  
  ggarrange(bar,clinical_TumorSize,#TIDE_bar3,
            ncol = 1,nrow = 2,align = "v",heights = c(100,10),legend = "none")
  
  ggarrange(clinical_Cancer,bar,TIDE_Exclusion,TIDE_CAF,clinical_Right,clinical_MMR,clinical_TumorSize,clinical_T,clinical_Stage,
            ncol = 1,nrow = 9,align = "v",heights = c(5,100,10,10,10,10,10,10,10),legend = "none")
  
  }
  
  }

#tissue
{
  Immune_fraction <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/res_cibersort_multiomics_plasma_PBMC_tissue_LM22_abs_clustered.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
  Immune_fraction <- Immune_fraction[,-grep("P-value|Correlation|RMSE|Absolute score",colnames(Immune_fraction))]
  #Immune_fraction <- Immune_fraction[-grep("CRC.PKU.15|STAD.PKU.23|STAD.PKU.36",rownames(Immune_fraction)),]
  #Immune_fraction <- Immune_fraction[grep("STAD",rownames(Immune_fraction)),]
  Immune_fraction <- Immune_fraction[-which(rowSums(Immune_fraction)==0),]
  Immune_fraction <- Immune_fraction/rowSums(Immune_fraction)
  Immune_fraction$ID <- rownames(Immune_fraction)
  Immune_fraction$ID <- gsub(".","-",fixed = TRUE,Immune_fraction$ID)
  
  rank1 <- Immune_fraction[Immune_fraction$`T cells gamma delta`>0.0,]
  rank2 <- Immune_fraction[Immune_fraction$`T cells gamma delta`<=0.0,]
  final_rank <- c(rank1[order(rank1$`T cells gamma delta`+rank1$`B cells`+rank1$`Plasma cells`,decreasing = TRUE),]$ID,
                  rank2[order(rank2$`B cells`+rank2$`Plasma cells`,decreasing = TRUE),]$ID)
  
  
  Immune_fraction_forplot <- melt(Immune_fraction,id.vars = "ID")
  Immune_fraction_forplot <- Immune_fraction_forplot[-grep("NC|Normal|Tumor|pico|mix",Immune_fraction_forplot$ID),]
  
  Immune_fraction_forplot$ID <- factor(Immune_fraction_forplot$ID,levels = final_rank)
  #Immune_fraction_forplot$ID <- factor(Immune_fraction_forplot$ID,levels = Immune_fraction[order((Immune_fraction$`CD8 T cells`-Immune_fraction$`T cells gamma delta`),decreasing = FALSE),]$ID)
  
  {
    clinical <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Sample_clinical_information.csv",header = TRUE,check.names = FALSE,row.names = 1)
    clinical$ID <- paste0(rownames(clinical),"-Tumor")
    #clinical <- clinical[-grep("NC|PBMC|Tumor|Normal|mix",clinical$ID),]
    rownames(clinical) <- clinical$ID
    clinical <- clinical[as.character(unique(Immune_fraction_forplot$ID)),]
    clinical$ID <- factor(clinical$ID,levels=levels(Immune_fraction_forplot$ID))
    #clinical$ID <- factor(clinical$ID,levels=c(clinical[order(clinical$`Tumor size`,decreasing = TRUE),]$ID))
    #clinical$ID <- factor(clinical$ID,levels=c(clinical[order(clinical$T,decreasing = TRUE),]$ID))
    #Immune_fraction_forplot$ID <- factor(Immune_fraction_forplot$ID,levels = levels(clinical$ID))
  }
  
  Immune_fraction_forplot$variable <- factor(Immune_fraction_forplot$variable, levels = rev(c("T cells gamma delta","B cells","Plasma cells","CD8 T cells","CD4 T cells",
                                                                                              "NK cells","Monocytes","Macrophages","Dendritic cells",
                                                                                              "Mast cells","Eosinophils","Neutrophils")))
  simpsons <- c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF",
                "#197EC0FF",alpha("purple",alpha=0.4),alpha("brown",alpha=0.8),"#71D0F5FF","#370335FF","dark green",
                "#C80813FF","#91331FFF","#1A9993FF","#FD8CC1FF")
  bar <- 
    ggplot(Immune_fraction_forplot,aes(x=ID,y=value,fill=variable))+
    geom_bar(stat = "identity",position = "stack")+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(10,10,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      axis.text.x = element_text(color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
      axis.text.y = element_text(face="bold", color="black", size=16, angle = 90,hjust=0.5),
      axis.title.x = element_text(face="bold", color="black", size=20),
      axis.title.y = element_text(face="bold",color="black", size=20))+
    #scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","","0.5","","1"),expand = c(0,0.01),limits = c(0,1))+
    scale_fill_manual(values= simpsons)
    #scale_fill_simpsons(alpha = 1)
  
  
  
  clinical_bar1 <-
    ggplot(clinical,aes(x=ID,y=0,fill= log2(clinical$`Tumor size`)))+
    geom_tile()+
    scale_fill_gradient2(low = "white",
                         mid = "grey",
                         high = "black",
                         midpoint = log2(5))+
    #scale_fill_manual(values=c("-"="white","early"="sky blue","late"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_bar2 <-
    ggplot(clinical,aes(x=ID,y=0,fill= clinical$PLT))+#`Tumor size`))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "orange")+
    #scale_fill_manual(values=c("-"="white","early"="sky blue","late"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_bar3 <-
    ggplot(clinical,aes(x=ID,y=0,fill= clinical$P53))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","-"="white","+"="light pink","++"="pink","+++"="dark red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_bar4 <-
    ggplot(clinical,aes(x=ID,y=0,fill= clinical$`MSI/MMR`))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","Normal"="white","Unstable; lack"="dark red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_bar5 <-
    ggplot(clinical,aes(x=ID,y=0,fill= clinical$HER2))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","-"="white","+"="light pink","++"="pink","+++"="dark red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_bar6 <-
    ggplot(clinical,aes(x=ID,y=0,fill= clinical$`Right half / left half`))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","L"="white","R"="red","L+R"="dark red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_bar7 <-
    ggplot(clinical,aes(x=ID,y=0,fill= clinical$Stage))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","1"="white","1A"="white","1B"="white","2A"="white","2B"="white","2C"="white","3B"="red","3C"="red","4"="dark red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  clinical_bar8 <-
    ggplot(clinical,aes(x=ID,y=0,fill= clinical$T))+#`Tumor size`))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "black")+
    scale_fill_manual(values=c("NA"="grey","Tis"="grey","1"="white","1a"="white","1b"="white","2"="white","3"="red","4a"="dark red","4b"="dark red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  
  TIDE <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Tissue_TIDE.txt",sep="\t",header = TRUE,check.names = FALSE,row.names = 1)
  TIDE$ID <- gsub("T","Tumor",rownames(TIDE))
  rownames(TIDE) <-  TIDE$ID
  TIDE <- TIDE[as.character(unique(Immune_fraction_forplot$ID)),]
  TIDE$ID <- factor(TIDE$ID,levels(Immune_fraction_forplot$ID))
  
  TIDE_bar1 <-
    ggplot(TIDE,aes(x=ID,y=0,fill= Responder))+
    geom_tile()+
    #scale_fill_gradient(low="white",high = "blue")+
    scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  TIDE_bar2 <-
    ggplot(TIDE,aes(x=ID,y=0,fill= IFNG))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "sky blue")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  TIDE_bar3 <-
    ggplot(TIDE,aes(x=ID,y=0,fill= MDSC+CAF))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "blue")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    #title("MDSC+CAF")+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  EPIC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Tissue_EPIC_cellFraction_TRef.txt",sep="\t",header = TRUE,check.names = FALSE,row.names = 1)
  EPIC$ID <- gsub("T","Tumor",rownames(EPIC))
  rownames(EPIC) <-  EPIC$ID
  EPIC <- EPIC[as.character(unique(Immune_fraction_forplot$ID)),]
  EPIC$ID <- factor(EPIC$ID,levels(Immune_fraction_forplot$ID))
  
  EPIC_bar1 <-
    ggplot(EPIC,aes(x=ID,y=0,fill= EPIC$CAFs+EPIC$Endothelial))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "blue")+
    #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-50,0,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  
  ggarrange(bar,clinical_bar2,clinical_bar3,clinical_bar4,clinical_bar5,clinical_bar6,clinical_bar7,clinical_bar8,TIDE_bar1,TIDE_bar2,TIDE_bar3,EPIC_bar1,clinical_bar1,
            ncol = 1,nrow = 13,align = "v",heights = c(100,5,5,5,5,5,5,5,5,5,5,5,0.0005),legend = "none")
}

#get CAF genes
{
  #method1:from EPIC, elife-26476-supp2-v2, CAF expression / all cell expression > 0.99 genes
  library(biomaRt)
  genes <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/GeneList_from_ImmPort.csv",header=TRUE,check.names = FALSE)
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  gene_names <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol","chromosome_name"),
                      filters = "hgnc_symbol",
                      values=as.character(genes$Gene.Symbol), mart= mart,useCache = FALSE)
  
  write.csv(gene_names[-grep("CHR",gene_names$chromosome_name),],"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/genes_ensembl_immport.csv")
  
  #method2:from TIDE, and the reference 83 Calon, A. et al. Dependency of colorectal cancer on a TGF-beta-driven program in stromal cells for metastasis initiation. Cancer cell 22, 571-584, doi:10.1016/j.ccr.2012.08.013 (2012).
  library(affy)
  library(limma)
  library(biomaRt)
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  
  path_to_cel <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/GSE39396_RAW/"
  datafiles <- ReadAffy(celfile.path = path_to_cel)
  eset <- rma(datafiles)
  View(eset)
  write.exprs(eset,file="/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/GSE39396.txt")
  
  Expression_matrix <- data.frame(exprs(eset))
  ph<-pData(eset) 
  design <- model.matrix(~factor(p_disease)) 
  colnames(design) <- c("case","control") 
  
  fit <- lmFit(eset, design) 
  fit <- eBayes(fit) 
  options(digits=2) 
  
  genes<- topTable(fit, coef=2, n=40, adjust="BH") 
  
  gene_names <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol","chromosome_name"),
                      filters = "affy_ht_hg_u133_plus_pm",
                      values=rownames(Expression_matrix), mart= mart,useCache = FALSE)
  
}

#gene level correlation
{
  Plasma_TPM <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/TPM/Plasma_PBMC_tissue_TPM_noduplicated_gene.txt",sep = "\t",check.names = FALSE)
  GSE173402_TPM <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/TPM/GSE174302_intron-spanning-available-samples_TPM.txt",sep = "\t",check.names = FALSE)
  colnames(GSE173402_TPM) <- gsub(".","-",fixed=TRUE,colnames(GSE173402_TPM))
  Plasma_TPM$ID <- rownames(Plasma_TPM)
  GSE173402_TPM$ID <- rownames(GSE173402_TPM)
  
  Plasma_TPM <- left_join(Plasma_TPM,GSE173402_TPM, by = c("ID"="ID"))
  Plasma_TPM[is.na(Plasma_TPM)] <- 0
  rownames(Plasma_TPM) <- Plasma_TPM$ID
  Plasma_TPM <- Plasma_TPM[,-which(colnames(Plasma_TPM)=="ID")]
  
  #LM22_genes <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/signature/LM22_ENSG.txt",sep = "\t")
  #Marker_genes <- Plasma_TPM[as.character(LM22_genes$Gene.symbol),as.character(unique(Immune_fraction_forplot$ID))]
  
  Interested_genes <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/genes_ensembl_Multiomics_STAD_potential_genes.csv")
  Marker_genes <- Plasma_TPM[unique(as.character(Interested_genes$ensembl_gene_id)),as.character(unique(Immune_fraction_forplot$ID))]
  Marker_genes <- na.omit(Marker_genes)
  
  forCorrlation <- cbind(clinical,t(Marker_genes))
  
  forCorrlation <- forCorrlation[grep("pico",rownames(forCorrlation)),]
  #forCorrlation <- forCorrlation[which(forCorrlation$Type=="Colorectal Cancer"),]
  
  forCorrlation <- forCorrlation[which(forCorrlation$Type=="Stomach Cancer"),]
  forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="Stage")]
  forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="ID")]
  forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="N")]
  forCorrlation <- forCorrlation[,-grep("Type|Patient ID|Gender|Age|Stage_simplified|T|M|Diagnosis|Tumor size|Vascular tumor thrombus|Neurological invasion",colnames(forCorrlation))]
  forCorrlation <- forCorrlation[,-grep("Tumor deposition/cancer nodules|Multiple primary cancer|Hgb|PLT|ALT|AST|ALB|PT|CA724|CA242|CEA|CA199|EGFR|HER2|CDNA2|TOP IIa|Right half / left half|MMR|MLH1|PMS2|MSH2|MSH6|MUC2|MUC5AC|EBER",colnames(forCorrlation))]
  
  correlation_method <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Multiomics_STAD_genes/Multiomics_STAD_genes"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Multiomics_STAD_genes/")
  enmuerate_correlation(forCorrlation,correlation_method,output_dir)
}

#Fit correlation
{
  
  Interested_genes <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/genes_ensembl_CAF_EPIC.csv")
  Marker_genes <- Plasma_TPM[unique(as.character(Interested_genes$ensembl_gene_id)),as.character(unique(Immune_fraction_forplot$ID))]
  Marker_genes <- na.omit(Marker_genes)
  forCorrlation <- cbind(clinical,t(Marker_genes))
  
  #forfit <- forCorrlation[grep("pico",rownames(forCorrlation)),]
  #forfit <- forfit[grep("STAD",rownames(forfit)),]
  #forvalidation <- forCorrlation[-grep("pico",rownames(forCorrlation)),]
  #forvalidation <- forvalidation[grep("STAD",rownames(forvalidation)),]
  forfit <- forCorrlation[grep("pico",rownames(forCorrlation)),]
  forfit <- forfit[grep("STAD|CRC",rownames(forfit)),]
  forvalidation <- forCorrlation[-grep("pico",rownames(forCorrlation)),]
  forvalidation <- forvalidation[grep("STAD|CRC",rownames(forvalidation)),]
  
  {
  #method1
  library(glmnet)
  #CAFs
  A <- as.matrix(sapply(forfit[,c("Stage_raw","ENSG00000108821","ENSG00000171345","ENSG00000168542","ENSG00000162366","ENSG00000164692","ENSG00000172426","ENSG00000141448","ENSG00000091986","ENSG00000174939","ENSG00000153495","ENSG00000115380","ENSG00000178531","ENSG00000186832","ENSG00000155622","ENSG00000006042","ENSG00000060718","ENSG00000175745","ENSG00000205076","ENSG00000183569","ENSG00000047936")], as.numeric))
  #gamma delta
  #A <- as.matrix(sapply(forfit[,c("Stage_raw","ENSG00000089012","ENSG00000183918","ENSG00000271503","ENSG00000211751","ENSG00000172116","ENSG00000125245","ENSG00000211829","ENSG00000174946","ENSG00000167286","ENSG00000206561","ENSG00000158050","ENSG00000139187","ENSG00000077984","ENSG00000089692","ENSG00000145649","ENSG00000112303","ENSG00000115607","ENSG00000113088")], as.numeric))
  #gamma delta
  #A <- as.matrix(sapply(forfit[,c("Stage_raw","ENSG00000089012","ENSG00000122224","ENSG00000162676","ENSG00000183918","ENSG00000213658")], as.numeric))
  #Multiomics_CRC_genes_correlated with stage, pvalue < 0.05
  A <- as.matrix(sapply(forfit[,c("Stage_raw","ENSG00000145088","ENSG00000019169","ENSG00000100095","ENSG00000172137",
                                  "ENSG00000168542","ENSG00000282608","ENSG00000135094","ENSG00000214140","ENSG00000175920",
                                  "ENSG00000140932","ENSG00000142208","ENSG00000100122","ENSG00000162676","ENSG00000082701",
                                  "ENSG00000230489","ENSG00000165025","ENSG00000104432","ENSG00000161911","ENSG00000008441",
                                  "ENSG00000166501","ENSG00000177885","ENSG00000112818","ENSG00000101096","ENSG00000204193",
                                  "ENSG00000104894","ENSG00000243156","ENSG00000174775","ENSG00000129194","ENSG00000145423",
                                  "ENSG00000122870","ENSG00000117020")], as.numeric))
  
  #Multiomics_STAD_genes_correlated with stage, pvalue < 0.05
  A <- as.matrix(sapply(forfit[,c("Stage_raw","ENSG00000163735","ENSG00000176547","ENSG00000140749","ENSG00000081059",
                                  "ENSG00000167286","ENSG00000213413","ENSG00000271503","ENSG00000221882","ENSG00000173372",
                                  "ENSG00000026559","ENSG00000282608","ENSG00000137462","ENSG00000164047","ENSG00000136689",
                                  "ENSG00000100721","ENSG00000123405","ENSG00000163563","ENSG00000197249","ENSG00000003137",
                                  "ENSG00000100362","ENSG00000110077","ENSG00000173083","ENSG00000163734","ENSG00000204472",
                                  "ENSG00000158473","ENSG00000167414","ENSG00000100450","ENSG00000188859","ENSG00000250510",
                                  "ENSG00000227507","ENSG00000158874","ENSG00000158825","ENSG00000166289","ENSG00000172548",
                                  "ENSG00000171051","ENSG00000213402","ENSG00000132514")], as.numeric))
  
  #fit_lasso <- glmnet(A[,-1] , A[,1] , standardize = TRUE , alpha =0.5) #LASSO model
  #print(fit_lasso) #LASSO model for different lambdas
  
  cvfit <- cv.glmnet( A[,-1] , A[,1] , standardize = TRUE , type.measure = "mse" , nfolds = 5 , alpha = 1) 
  cvfit    
  cvfit$lambda.min
  cvfit$lambda.1se
  coef(cvfit , s = "lambda.min")
  coef(cvfit , s = "lambda.1se")
  
  #method2 ?? not certain
  {
  library(tidyverse)
  library(caret)
  A <- as.matrix(sapply(forCorrlation[,c("Stage_raw","ENSG00000108821","ENSG00000171345","ENSG00000168542","ENSG00000162366","ENSG00000164692","ENSG00000172426","ENSG00000141448","ENSG00000091986","ENSG00000174939","ENSG00000153495","ENSG00000115380","ENSG00000178531","ENSG00000186832","ENSG00000155622","ENSG00000006042","ENSG00000060718","ENSG00000175745","ENSG00000205076","ENSG00000183569","ENSG00000047936")], as.numeric))
  A$`Tumor size` <- as.numeric(as.character(A$`Tumor size`))
  A$`T cells gamma delta` <- as.numeric(as.character(A$`T cells gamma delta`))
  A$`T cells CD4` <- as.numeric(as.character(A$`T cells CD4`))
  A$Exclusion <- as.numeric(as.character(A$Exclusion))
  A$CAF <- as.numeric(as.character(A$CAF))
  A$`B cells` <- as.numeric(as.character(A$`B cells`))
  A$`Plasma cells` <- as.numeric(as.character(A$`Plasma cells`))
  head(A)
  
  set.seed(1998) # 设置随机种子使结果可重复
  idx = sample(nrow(A), 0.8 * nrow(A))
  trainData  <- A[idx, ]
  testData   <- A[-idx, ]
  
  model <- lm(`Tumor size` ~ `T cells gamma delta` + `T cells CD4` + `Exclusion` + `CAF`, 
              data = trainData)
  
  summary(model)$coef
  
  predict(model,testData)
  
  model <- lm(`Tumor size` ~ `Exclusion` + `T cells gamma delta` + `B cells` + `Plasma cells`, 
              data = A)
  
  summary(model)$coef
  }
}

#CAFs
forfit$`Fitted` <- 0.2084667051*forfit$ENSG00000108821+0.0699843508*forfit$ENSG00000171345+0.3286087154*forfit$ENSG00000168542-0.0002842341*forfit$ENSG00000162366+ 
0.4808993020*forfit$ENSG00000164692-0.0673166255*forfit$ENSG00000172426+0.0889751549*forfit$ENSG00000141448+0.2194534345*forfit$ENSG00000091986+
0.6918664090*forfit$ENSG00000174939+0.0339712745*forfit$ENSG00000153495+0.0429682274*forfit$ENSG00000115380-0.5139937159*forfit$ENSG00000178531+ 
0.0243741768*forfit$ENSG00000186832-0.0472012043*forfit$ENSG00000155622+0.1070337255*forfit$ENSG00000006042+0.1790269326*forfit$ENSG00000060718+  
0.2421727418*forfit$ENSG00000175745+0.5693757555*forfit$ENSG00000205076-0.1713988864*forfit$ENSG00000183569+0.3891344341*forfit$ENSG00000047936  

forfit$`Fitted` <- 
  0.209437207*forfit$ENSG00000171345+0.742854743*forfit$ENSG00000164692-0.009181138*forfit$ENSG00000172426+0.074091423*forfit$ENSG00000115380-0.061449374*forfit$ENSG00000155622 
#gammadeltaT 
forfit$`Fitted` <- -6.059254e-04*forfit$ENSG00000089012-5.665177e-04*forfit$ENSG00000183918-1.280746e-06*forfit$ENSG00000271503-8.250298e-05*forfit$ENSG00000211751+
1.162115e-04*forfit$ENSG00000172116-6.216530e-04*forfit$ENSG00000125245-2.667684e-05*forfit$ENSG00000211829-2.455918e-03*forfit$ENSG00000174946 
-7.147134e-05*forfit$ENSG00000167286-7.103221e-04*forfit$ENSG00000206561-3.057003e-03*forfit$ENSG00000158050-1.694696e-04*forfit$ENSG00000139187 
-2.824135e-05*forfit$ENSG00000077984-7.939934e-04*forfit$ENSG00000089692-4.949101e-05*forfit$ENSG00000145649-1.090442e-04*forfit$ENSG00000112303 
-1.743048e-04*forfit$ENSG00000115607+7.649479e-05*forfit$ENSG00000113088


#Multiomics_CRC_plasma_score 
#check on train set (multiomics data, 5-fold, alpha = 0.5)
forfit$`Fitted` <- -0.013119497*forfit$ENSG00000145088+0.032007408*forfit$ENSG00000100095-4.119734563*forfit$ENSG00000282608-2.887588405*forfit$ENSG00000135094+3.101657587*forfit$ENSG00000214140-1.406268603*forfit$ENSG00000175920+0.006986278*forfit$ENSG00000140932-0.372400199*forfit$ENSG00000100122-0.157786357*forfit$ENSG00000230489+0.001228865*forfit$ENSG00000161911+2.809393881*forfit$ENSG00000112818+0.183997173*forfit$ENSG00000174775+2.673318410*forfit$ENSG00000129194-0.865182732*forfit$ENSG00000122870-0.005341222*forfit$ENSG00000117020 

#check on train set (multiomics data, 5-fold, alpha = 0.9)
forfit$`Fitted` <- -0.014383048*forfit$ENSG00000145088-4.139422564*forfit$ENSG00000282608-2.587594386*forfit$ENSG00000135094+2.241242197*forfit$ENSG00000214140-1.154165765*forfit$ENSG00000175920+0.004881598*forfit$ENSG00000140932-0.184916236*forfit$ENSG00000100122-0.133444656*forfit$ENSG00000230489+0.001148320*forfit$ENSG00000161911+2.610687832*forfit$ENSG00000112818+0.160862183*forfit$ENSG00000174775+2.323295412*forfit$ENSG00000129194-0.391285458*forfit$ENSG00000122870-0.003163634*forfit$ENSG00000117020 

#Multiomics_STAD_plasma_score
#check on train set (multiomics data, 5-fold, alpha = 1)
forfit$`Fitted` <- 0.002645357*forfit$ENSG00000176547-0.036576575*forfit$ENSG00000164047 


forfit$Stage_raw <- as.numeric(forfit$Stage_raw)
ExclusionScore_TumorSize <-
  ggscatter(forfit, x = "Fitted", y = "Stage_raw", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
  #xlab("Exclusion - γδ T cells - B cells - plasma cells")+
  #ylab("Tumor size (cm)")+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    plot.margin = unit(c(20,20,20,20),"pt"),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=20, family = "Arial"),
    axis.text.y = element_text(face="bold",  color="black", size=20, family = "Arial"),
    axis.title.x = element_text(face="bold", color="black", size=16, family = "Arial"),
    axis.title.y = element_text(face="bold",color="black", size=24, family = "Arial"))
ExclusionScore_TumorSize
ggsave(plot = ExclusionScore_TumorSize, filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Multiomics_fitted_STAD_plasma_score.pdf",device = "pdf",width = 5.75,height = 4.65)

#check on external data (GSE174302)
#alpha = 0.5
forvalidation$`Fitted` <- -0.013119497*forvalidation$ENSG00000145088+0.032007408*forvalidation$ENSG00000100095-4.119734563*forvalidation$ENSG00000282608-2.887588405*forvalidation$ENSG00000135094+3.101657587*forvalidation$ENSG00000214140-1.406268603*forvalidation$ENSG00000175920+0.006986278*forvalidation$ENSG00000140932-0.372400199*forvalidation$ENSG00000100122-0.157786357*forvalidation$ENSG00000230489+0.001228865*forvalidation$ENSG00000161911+2.809393881*forvalidation$ENSG00000112818+0.183997173*forvalidation$ENSG00000174775+2.673318410*forvalidation$ENSG00000129194-0.865182732*forvalidation$ENSG00000122870-0.005341222*forvalidation$ENSG00000117020 
#alpha = 0.9
forvalidation$`Fitted` <- -0.014383048*forvalidation$ENSG00000145088-4.139422564*forvalidation$ENSG00000282608-2.587594386*forvalidation$ENSG00000135094+2.241242197*forvalidation$ENSG00000214140-1.154165765*forvalidation$ENSG00000175920+0.004881598*forvalidation$ENSG00000140932-0.184916236*forvalidation$ENSG00000100122-0.133444656*forvalidation$ENSG00000230489+0.001148320*forvalidation$ENSG00000161911+2.610687832*forvalidation$ENSG00000112818+0.160862183*forvalidation$ENSG00000174775+2.323295412*forvalidation$ENSG00000129194-0.391285458*forvalidation$ENSG00000122870-0.003163634*forvalidation$ENSG00000117020 

#Multiomics_STAD_plasma_score
#check on external set (GSE174302)
forvalidation$`Fitted` <- 0.002645357*forvalidation$ENSG00000176547-0.036576575*forvalidation$ENSG00000164047 

forvalidation$`Fitted` <- 
  0.209437207*forvalidation$ENSG00000171345+0.742854743*forvalidation$ENSG00000164692-0.009181138*forvalidation$ENSG00000172426+0.074091423*forvalidation$ENSG00000115380-0.061449374*forvalidation$ENSG00000155622 

forvalidation$Stage_raw <- as.numeric(forvalidation$Stage_raw)
ExclusionScore_TumorSize <-
  ggscatter(forvalidation, x = "Fitted", y = "Stage_raw", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
  #xlab("Exclusion - γδ T cells - B cells - plasma cells")+
  #ylab("Tumor size (cm)")+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    plot.margin = unit(c(20,20,20,20),"pt"),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=20, family = "Arial"),
    axis.text.y = element_text(face="bold",  color="black", size=20, family = "Arial"),
    axis.title.x = element_text(face="bold", color="black", size=16, family = "Arial"),
    axis.title.y = element_text(face="bold",color="black", size=24, family = "Arial"))
ExclusionScore_TumorSize
ggsave(plot = ExclusionScore_TumorSize, filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Multiomics_validationonGSE174302_CRC_plasma_score.pdf",device = "pdf",width = 5.75,height = 4.65)

fitted_order <- rownames(forCorrlation[order(forCorrlation$Exclusion,decreasing = TRUE),])

#Immune_fraction_forplot$ID <- factor(Immune_fraction_forplot$ID,levels = fitted_order)
#clinical$ID <- factor(clinical$ID,levels=fitted_order)
#TIDE$ID <- factor(TIDE$ID,levels=fitted_order)
#EPIC$ID <- factor(EPIC$ID,levels = fitted_order)
}

#Survival analysis & stage analysis
{
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/21.Survival analysis")
  library(survival)
  library(survminer)
  library(TCGAutils)
  library(dplyr)
  
  #read sample sheet
  sample_sheet1 <- read.table('COAD_TCGA/COAD_sample_sheet.2019-06-08.tsv',header=T,sep='\t',check.names = FALSE)
  sample_sheet2 <- read.table('READ_TCGA/READ_sample_sheet.2019-06-08.tsv',header=T,sep='\t',check.names = FALSE)
  sample_sheet3 <- read.table('STAD_TCGA/STAD_sample_sheet.2019-06-09.tsv',header=T,sep='\t',check.names = FALSE)
  sample_sheet <- rbind(sample_sheet1,sample_sheet2,sample_sheet3)
  sample_sheet <- sample_sheet[grep("htseq.counts.gz",sample_sheet$`File Name`),]
  sample_sheet$`File Name` <- as.character(lapply(strsplit(as.character(sample_sheet$`File Name`),".",fixed = TRUE), function(x) x[1]))
  
  head(sample_sheet)
  
  # read RNA file 
  rna1 <- read.table('COAD_TCGA/TCGA_COAD_count_matrix_TPM.txt',header=T,row.names=1,sep='\t',check.names = FALSE)
  rna2 <- read.table('READ_TCGA/TCGA_READ_count_matrix_TPM.txt',header=T,row.names=1,sep='\t',check.names = FALSE)
  rna3 <- read.table('STAD_TCGA/TCGA_STAD_count_matrix_TPM.txt',header=T,row.names=1,sep='\t',check.names = FALSE)
  rna <- cbind(rna1,rna2,rna3)
  rna_filename <- data.frame("file_name"=colnames(rna))
  rna_caseid <- left_join(rna_filename,sample_sheet,by=c("file_name"="File Name"))
  
  # and read the Clinical file, in this case i transposed it to keep the clinical feature title as column name
  clinical1 <- t(read.table('COAD_TCGA/clinical.tsv',header=T, row.names=1, sep='\t'))
  clinical2 <- t(read.table('READ_TCGA/clinical.tsv',header=T, row.names=1, sep='\t'))
  clinical3 <- t(read.table('STAD_TCGA/clinical.tsv',header=T, row.names=1, sep='\t'))
  clinical <- cbind(clinical1,clinical2,clinical3)
  clinical_caseid <- UUIDtoBarcode(colnames(clinical), from_type = "case_id")
  
  complete_id <- left_join(rna_caseid,clinical_caseid,by=c("Case ID"="submitter_id"))
  complete_id <- complete_id[complete_id$`Sample Type`=="Primary Tumor",]
  
  #substitute
  rna_forSurvival <- rna[,complete_id[grep("TCGA-STAD|TCGA-READ|TCGA-COAD",complete_id$`Project ID`),]$file_name]
  colnames(rna_forSurvival) <- complete_id[grep("TCGA-STAD|TCGA-READ|TCGA-COAD",complete_id$`Project ID`),]$`Sample ID`
  clinical_forSurvival <- clinical[,complete_id[grep("TCGA-STAD|TCGA-READ|TCGA-COAD",complete_id$`Project ID`),]$case_id]
  colnames(clinical_forSurvival) <- complete_id[grep("TCGA-STAD|TCGA-READ|TCGA-COAD",complete_id$`Project ID`),]$`Sample ID`
  
  clinical_forSurvival.t <- as.data.frame(t(clinical_forSurvival))
  No_status <- paste(names(grep("--|Not Reported",clinical_forSurvival.t$vital_status,value = TRUE)),collapse = "|")
  complete_id <- complete_id[-grep(No_status,complete_id$`Sample ID`),]
  #rm(rna)
  #rm(clinical)
  
  # first I remove genes whose expression is == 0 in more than 50% of the samples:
  rem <- function(x){
    x <- as.matrix(x)
    x <- t(apply(x,1,as.numeric))
    r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
    remove <- which(r > dim(x)[2]*0.5)
    return(remove)
  }
  remove <- rem(rna_forSurvival)
  rna_forSurvival <- rna_forSurvival[-remove,]
  
  
  #Now I need to identify normal and tumor samples. this is done using the TCGA barcode (https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode). The two digits at position 14-15 of the barcode will indicate teh sample type, from the link:
  #  "Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29."
  
  # see the values
  table(substr(colnames(rna_forSurvival),14,14))
  # 0      1 
  # 534   72
  
  clinical_forSurvival.t <- as.data.frame(t(clinical_forSurvival))
  #for alive indvivduals, time = days to last follow up; for dead indvivuals, time = days to death
  clinical_forSurvival.t$time <- NA
  clinical_forSurvival.t[clinical_forSurvival.t$vital_status=="Alive",]$time <- as.numeric(as.character(clinical_forSurvival.t[clinical_forSurvival.t$vital_status=="Alive",]$days_to_last_follow_up))
  clinical_forSurvival.t[clinical_forSurvival.t$vital_status=="Dead",]$time <- as.numeric(as.character(clinical_forSurvival.t[clinical_forSurvival.t$vital_status=="Dead",]$days_to_death))
  
  rna_forSurvival.t <- as.data.frame(t(rna_forSurvival))
  
  forSurvival <- cbind(clinical_forSurvival.t,rna_forSurvival.t)
  forSurvival$days_to_death <- as.numeric(as.character(forSurvival$days_to_death))
  forSurvival$vital_status <- as.character(forSurvival$vital_status)
  forSurvival$vital_status <- gsub("Alive","1",forSurvival$vital_status)
  forSurvival$vital_status <- gsub("Dead","2",forSurvival$vital_status)
  forSurvival$vital_status <- as.integer(forSurvival$vital_status)
  
  forSurvival$Customized_group <- NA
  
  #grep("TP53",colnames(forSurvival),value = TRUE)
  #forSurvival[forSurvival[["ENSG00000141510|TP53|2653"]] >= median(forSurvival[["ENSG00000141510|TP53|2653"]]),]$Customized_group <- paste0("ENSG00000141510|TP53|2653", "|High expressed group (50%)")
  #forSurvival[forSurvival[["ENSG00000141510|TP53|2653"]] < median(forSurvival[["ENSG00000141510|TP53|2653"]]),]$Customized_group <- paste0("ENSG00000141510|TP53|2653", "|Low expressed group (50%)")
  
  #CAFs
  forSurvival$Plasma_exclusion_score <-  0.2084667051*forSurvival$`ENSG00000108821|COL1A1|5914`+0.0699843508*forSurvival$`ENSG00000171345|KRT19|1390`+0.3286087154*forSurvival$`ENSG00000168542|COL3A1|5490`-0.0002842341*forSurvival$`ENSG00000162366|PDZK1IP1|838`+ 
    0.4808993020*forSurvival$`ENSG00000164692|COL1A2|5993`-0.0673166255*forSurvival$`ENSG00000172426|RSPH9|2586`+0.0889751549*forSurvival$`ENSG00000141448|GATA6|3624`+0.2194534345*forSurvival$`ENSG00000091986|CCDC80|12301`+
    0.6918664090*forSurvival$`ENSG00000174939|ASPHD1|1816`+0.0339712745*forSurvival$`ENSG00000153495|TEX29|835`+0.0429682274*forSurvival$`ENSG00000115380|EFEMP1|3024`-0.5139937159*forSurvival$`ENSG00000178531|CTXN1|1237`+ 
    0.0243741768*forSurvival$`ENSG00000186832|KRT16|1658`-0.0472012043*forSurvival$`ENSG00000155622|XAGE2|627`+0.1070337255*forSurvival$`ENSG00000006042|TMEM98|4218`+0.1790269326*forSurvival$`ENSG00000060718|COL11A1|7327`+  
    0.2421727418*forSurvival$`ENSG00000175745|NR2F1|3843`+0.5693757555*forSurvival$`ENSG00000205076|LGALS7|565`-0.1713988864*forSurvival$`ENSG00000183569|SERHL2|2096`+0.3891344341*forSurvival$`ENSG00000047936|ROS1|8451`
  
  forSurvival$Plasma_exclusion_score <-  forSurvival$`ENSG00000108821|COL1A1|5914`+forSurvival$`ENSG00000171345|KRT19|1390`+forSurvival$`ENSG00000168542|COL3A1|5490`+forSurvival$`ENSG00000162366|PDZK1IP1|838`+ 
    forSurvival$`ENSG00000164692|COL1A2|5993`+forSurvival$`ENSG00000172426|RSPH9|2586`+forSurvival$`ENSG00000141448|GATA6|3624`+ forSurvival$`ENSG00000091986|CCDC80|12301`+
    forSurvival$`ENSG00000174939|ASPHD1|1816`+forSurvival$`ENSG00000153495|TEX29|835`+forSurvival$`ENSG00000115380|EFEMP1|3024`+forSurvival$`ENSG00000178531|CTXN1|1237`+ 
    forSurvival$`ENSG00000186832|KRT16|1658`+forSurvival$`ENSG00000155622|XAGE2|627`+forSurvival$`ENSG00000006042|TMEM98|4218`+forSurvival$`ENSG00000060718|COL11A1|7327`+  
    forSurvival$`ENSG00000175745|NR2F1|3843`+forSurvival$`ENSG00000205076|LGALS7|565`+forSurvival$`ENSG00000183569|SERHL2|2096`+forSurvival$`ENSG00000047936|ROS1|8451`
  
  forSurvival$Plasma_exclusion_score <- forSurvival$`ENSG00000108821|COL1A1|5914`
  forSurvival$Plasma_exclusion_score <- forSurvival$`ENSG00000153563|CD8A|3048`
  forSurvival[forSurvival[["Plasma_exclusion_score"]] >= median(forSurvival[["Plasma_exclusion_score"]]),]$Customized_group <- "High" #paste0("Plasma_exclusion_score", "|High expressed group (50%)")
  forSurvival[forSurvival[["Plasma_exclusion_score"]] < median(forSurvival[["Plasma_exclusion_score"]]),]$Customized_group <- "Low" #paste0("Plasma_exclusion_score", "|Low expressed group (50%)")
  
  #gammadeltaT 
  forSurvival$Plasma_gammadeltaT_score <- -6.059254e-04*forSurvival$`ENSG00000089012|SIRPG|1732`-5.665177e-04*forSurvival$`ENSG00000183918|SH2D1A|2554`-1.280746e-06*forSurvival$`ENSG00000271503|CCL5|1352`-8.250298e-05*forSurvival$`ENSG00000211751|TRBC1|760`+
    1.162115e-04*forSurvival$`ENSG00000172116|CD8B|4794`-6.216530e-04*forSurvival$`ENSG00000125245|GPR18|1875`-2.667684e-05*forSurvival$`ENSG00000211829|TRDC|720`-2.455918e-03*forSurvival$`ENSG00000174946|GPR171|1810`-
    7.147134e-05*forSurvival$`ENSG00000167286|CD3D|701`-7.103221e-04*forSurvival$`ENSG00000206561|COLQ|6346`-3.057003e-03*forSurvival$`ENSG00000158050|DUSP2|1685`-1.694696e-04*forSurvival$`ENSG00000139187|KLRG1|1874`-
    2.824135e-05*forSurvival$`ENSG00000077984|CST7|891`-7.939934e-04*forSurvival$`ENSG00000089692|LAG3|2587`-4.949101e-05*forSurvival$`ENSG00000145649|GZMA|896`-1.090442e-04*forSurvival$`ENSG00000112303|VNN2|2014`-
    1.743048e-04*forSurvival$`ENSG00000115607|IL18RAP|2773`+7.649479e-05*forSurvival$`ENSG00000113088|GZMK|1506`
  
  #gammadeltaT2
  forSurvival$Plasma_gammadeltaT_score <- -0.053551401*forSurvival$`ENSG00000089012|SIRPG|1732`+0.170103962*forSurvival$`ENSG00000122224|LY9|5083`+0.140654566*forSurvival$`ENSG00000162676|GFI1|4554`-0.042064351*forSurvival$`ENSG00000183918|SH2D1A|2554`+0.002171451*forSurvival$`ENSG00000213658|LAT|2443`
  
  forSurvival$Plasma_gammadeltaT_score <- -forSurvival$`ENSG00000089012|SIRPG|1732`+forSurvival$`ENSG00000122224|LY9|5083`+forSurvival$`ENSG00000162676|GFI1|4554`-forSurvival$`ENSG00000183918|SH2D1A|2554`+forSurvival$`ENSG00000213658|LAT|2443`
  
  #forSurvival$Plasma_gammadeltaT_score <- forSurvival$`ENSG00000172116|CD8B|4794`
  forSurvival[forSurvival[["Plasma_gammadeltaT_score"]] >= median(forSurvival[["Plasma_gammadeltaT_score"]]),]$Customized_group <- "High" #paste0("Plasma_gammadeltaT_score", "|High expressed group (50%)")
  forSurvival[forSurvival[["Plasma_gammadeltaT_score"]] < median(forSurvival[["Plasma_gammadeltaT_score"]]),]$Customized_group <- "Low" #paste0("Plasma_gammadeltaT_score", "|Low expressed group (50%)")
  forSurvival$Customized_group <- factor(forSurvival$Customized_group, levels = c("Low","High"))
  #forSurvival <- forSurvival[-which(forSurvival$project_id=="TCGA-STAD"),]
  
  fit1 <- surv_fit(Surv(time, vital_status) ~ Customized_group, data = forSurvival[,c("gender","time","vital_status","project_id","Customized_group")])
  plot1 <- ggsurvplot(fit1, data = forSurvival, pval = TRUE, pval.method = TRUE, conf.int = TRUE, 
                      #risk.table = "absolute", 
                      palette = "aaas", 
                      ggtheme = theme_survminer(base_size = 16, base_family = "Arial",
                                                font.main = c(40,"plain","black"), font.tickslab = c(16,"plain","black"),
                                                font.x = c(24,"plain","black"),font.y = c(24,"plain","black"),
                                                font.legend = c(16,"plain","black"), legend = "right"))
  plot1
  ggsave(plot = plot1 , filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/21.Survival analysis/Plasma_gammadeltaT_score_COAD+READ_641.pdf",device = "pdf",width = 4.5,height = 5)
  
  ggscatter(forSurvival, x = "tumor_stage", y = "Plasma_exclusion_score", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
    #xlab("Exclusion - γδ T cells - B cells - plasma cells")+
    #ylab("Tumor size (cm)")+
    theme(#legend.position="right",
      legend.position="right",
      panel.grid=element_blank(),
      plot.margin = unit(c(20,20,20,20),"pt"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=20, family = "Arial"),
      axis.text.y = element_text(face="bold",  color="black", size=20, family = "Arial"),
      axis.title.x = element_text(face="bold", color="black", size=16, family = "Arial"),
      axis.title.y = element_text(face="bold",color="black", size=24, family = "Arial"))
  
  forSurvival$Stage_simplified <- NA
  forSurvival[which(forSurvival$tumor_stage=="not reported"),]$Stage_simplified <- "Not reported"
  forSurvival[which(forSurvival$tumor_stage=="--"),]$Stage_simplified <- "Not reported"
  
  forSurvival[which(forSurvival$tumor_stage=="stage i"),]$Stage_simplified <- "Stage I"
  forSurvival[which(forSurvival$tumor_stage=="stage ia"),]$Stage_simplified <- "Stage I"
  forSurvival[which(forSurvival$tumor_stage=="stage ib"),]$Stage_simplified <- "Stage I"
  
  forSurvival[which(forSurvival$tumor_stage=="stage ii"),]$Stage_simplified <- "Stage II"
  forSurvival[which(forSurvival$tumor_stage=="stage iia"),]$Stage_simplified <- "Stage II"
  forSurvival[which(forSurvival$tumor_stage=="stage iib"),]$Stage_simplified <- "Stage II"
  forSurvival[which(forSurvival$tumor_stage=="stage iic"),]$Stage_simplified <- "Stage II"
  
  forSurvival[which(forSurvival$tumor_stage=="stage iii"),]$Stage_simplified <- "Stage III"
  forSurvival[which(forSurvival$tumor_stage=="stage iiia"),]$Stage_simplified <- "Stage III"
  forSurvival[which(forSurvival$tumor_stage=="stage iiib"),]$Stage_simplified <- "Stage III"
  forSurvival[which(forSurvival$tumor_stage=="stage iiic"),]$Stage_simplified <- "Stage III"
  
  forSurvival[which(forSurvival$tumor_stage=="stage iv"),]$Stage_simplified <- "Stage IV"
  forSurvival[which(forSurvival$tumor_stage=="stage iva"),]$Stage_simplified <- "Stage IV"
  forSurvival[which(forSurvival$tumor_stage=="stage ivb"),]$Stage_simplified <- "Stage IV"
  
  forSurvival <- forSurvival[which(forSurvival$Stage_simplified!="Not reported"),]
  p<-ggplot(forSurvival,aes(x=forSurvival$Stage_simplified,y=log10(forSurvival$Plasma_gammadeltaT_score),fill = forSurvival$Stage_simplified))+
    geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
    geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 0.8),color = alpha("black",alpha = 0.2))+
    #scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","early"="dark green","late"="black"))+
    scale_fill_manual(values=c("Stage I"="#90E0EF","Stage II"="#00B4D8","Stage III"="#03045E","Stage III/IV"="#03045E","Stage IV"="black"))+
    #scale_fill_brewer(palette="Blues") +
    #ylim(0,25)+
    theme_bw()+
    #ylab(colnames(plot)[1])+
    theme(#legend.position="right",
      legend.position="none",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(size=1, colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      #axis.text.x = element_blank(),
      axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1),
      axis.text.y = element_text(color="black", size=24),
      axis.title.x = element_text(color="black", size=0),
      axis.title.y = element_text(color="black", size=16))+
    stat_compare_means(comparisons = 
                         list(c("Stage I","Stage IV"),c("Stage I","Stage III"),c("Stage I","Stage II")),
                       method = "wilcox.test",
                       method.args = list(alternative = "greater",paired = TRUE),
                       label = "p.signif",
                       size = 10,
                       vjust = 0.5)
  ggsave(plot=p,filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/Gammadelta_inTCGA.pdf",width=4.40,height = 5.69)
q}

############################################################################################
#Logistic regression to evaluate different gene signatures in early stage vs. late stage (binomial clinical data)
{
  #variable***
  gene_signature_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Gene_signatures"
  workdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Fitted_signature_forStage"
  dir.create(workdir)
  gene_signatures <- gsub(".csv","",dir(gene_signature_dir))
  #gene_signatures <- gene_signatures[-grep("Bacterial invasion of epithelial cells_KEGG",gene_signatures)]
  #gene_signatures <- gene_signatures[grep("MolecularCancer",gene_signatures)]
  #gene_signatures <- c("genes_ensembl_gammadeltT_LM22","genes_ensembl_CAF_EPIC")
  
  n=1
  result <- as.data.frame(matrix(numeric(0),ncol=6))
  colnames(result) <- c("Genes","Formula","AUROC of Training fitted score","AUROC of Training summed score","AUROC of Testing Fitted score in GSE174302","AUROC of Testing summed score in GSE174302")
  while(n<=length(gene_signatures)){
  message(n)
  dir.create(paste0(workdir,"/",gene_signatures[n]))
  setwd(paste0(workdir,"/",gene_signatures[n]))
  
  Interested_genes <- read.csv(paste0(gene_signature_dir,"/",gene_signatures[n],".csv"))
  clinical_data <- "Stage_simplified"
  negative <- "early"
  positive <- "late"
  
  #read in genes of different signatures for training
  {
    TPM_forfit <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/TPM/Plasma_PBMC_tissue_TPM_noduplicated_gene.txt",header = TRUE,check.names = FALSE,row.names = 1,sep = "\t") 
    TPM_forfit <- TPM_forfit[,-grep("NC|PBMC|-T|-N|HCC|LUAD|ESCA|Tumor|Normal|NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",colnames(TPM_forfit))]
    
    clinical <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Sample_clinical_information_x.csv",header = TRUE,check.names = FALSE,row.names = 1)
    rownames(clinical) <- paste0(rownames(clinical),"-pico")
    clinical <- clinical[-grep("NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",rownames(clinical)),]
    clinical <- clinical[colnames(TPM_forfit),]
    
    Marker_genes <- TPM_forfit[unique(as.character(Interested_genes$ensembl_gene_id)),]
    rownames(Marker_genes) <- unique(as.character(Interested_genes$ensembl_gene_id))
    Marker_genes[is.na(Marker_genes)] <- 0 
    #Marker_genes <- Marker_genes[-grep("NA",rownames(Marker_genes)),]
    Expression_clinical <- cbind(clinical,t(Marker_genes))
    Expression_clinical <- Expression_clinical[-grep("x",Expression_clinical[[clinical_data]]),]
    
    forfit <- Expression_clinical
  }
  
  #fit a signature score based on multiomics RNA data, by cv.glmnet
  {
  #A is a subset of forfit, A only contains interested 1 clinical response and multiple genes
  A <- as.matrix(sapply(forfit[,c(clinical_data,unique(as.character(Interested_genes$ensembl_gene_id)))], as.numeric))
  rownames(A) <- rownames(forfit)
    
  cvfit <- cv.glmnet(A[,-1], A[,1] , standardize = TRUE , type.measure = "auc", family = "binomial" , nfolds = 5 , alpha = 1) 
  cvfit    
  cvfit$lambda.min
  cvfit$lambda.1se
  coefficient <- coef(cvfit , s = "lambda.min")
  #coefficient <- coef(cvfit , s = "lambda.1se")
  
  if(length((coefficient@i+1)[-1])==0){
    message("Skipped: ",gene_signatures[n],"\nBecause coefficient number = 0")
    n=n+1
    next
    }
  Lasso_selected_genes <- as.data.frame(A[,(coefficient@i+1)[-1]])
  genes_coefficients <- coefficient@x[-1]
  
  i=1
  Fitted_signature_forfit = 0
  formulas = paste0(clinical_data," ~ ")
  while(i<=ncol(Lasso_selected_genes)){
    Fitted_signature_forfit = Fitted_signature_forfit + genes_coefficients[i]*Lasso_selected_genes[i]
    formulas <- paste0(formulas,"+",genes_coefficients[i],"*",colnames(Lasso_selected_genes[i]),"")
    i=i+1
  }
  Fitted_signature_forfit <- sigmoid(Fitted_signature_forfit)
  colnames(Fitted_signature_forfit) <- "Fitted_signature"
  }
  
  
  #train summed signature and fitted score seperately by caret
  {
    library(caret)
    library(pROC)
    #Fitted score AUC in train set by tree model
    {
    train_labels <- forfit[[clinical_data]]
    train_labels <- factor(train_labels, levels = c(negative,positive))
    set.seed(100000)
    #fix model parameters, such as mtry, ntree
    tune_mtry <- tuneRF(Fitted_signature_forfit,train_labels) #tune mtry
    tune_mtry <- as.data.frame(tune_mtry)
    mtry_best <- tune_mtry[order(tune_mtry$OOBError,decreasing=FALSE),][1,1]
    
    #train RF model
    RF_model_fitted <- randomForest(Fitted_signature_forfit,train_labels, mtry = mtry_best)
    
    train_predicted <- as.data.frame(RF_model_fitted$votes)
    train_roc.curve <- roc(train_labels,train_predicted[,which(colnames(train_predicted)==positive)],levels=c(negative,positive),direction=c("<"))
    Fitted_signature_auc_in_trainset <- train_roc.curve$auc
    }
    
    #Fitted score AUC in train set by LR model
    #There must be two or more independent variables, or predictors, for a logistic regression.
    {
      #params.grid <- expand.grid(alpha = c(0,0.5,1),lambda = c(0,0.01,0.1,1))
      # alpha: relative weighting of L1 and L2 regularization
      # lambda: degree of regularization, see glmnet documentation for detail
      
      # train是caret封装的一个用于调参的函数
      # 和rfeControl类似，trainControl也返回一个列表，用于控制train函数的行为
      # method同样可以是"boot","cv"等，number参数控制了交叉验证的折数和bootstraping的重复次数
      #tr.ctrl <- trainControl(method="cv",number = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
      #cv.fitted <- train(Fitted_signature_forfit,train_labels,method="glmnet",family="binomial",metric = "ROC",tuneGrid = params.grid,preProcess = NULL,trControl = tr.ctrl )
    }
    
    
    #get summed gene set signature
    {
      Summed_signature <- as.data.frame(rowSums(log2(forfit[,rownames(Marker_genes)]+1)))
      colnames(Summed_signature) <- "Summed_signature"
    }
    #Summed score AUC in train set
    {
      train_labels <- forfit[[clinical_data]]
      train_labels <- factor(train_labels, levels = c(negative,positive))
      set.seed(100000)
      #fix model parameters, such as mtry, ntree
      tune_mtry <- tuneRF(Summed_signature,train_labels) #tune mtry
      tune_mtry <- as.data.frame(tune_mtry)
      mtry_best <- tune_mtry[order(tune_mtry$OOBError,decreasing=FALSE),][1,1]
      
      #train RF model use n_sample-1
      RF_model_summed <- randomForest(Summed_signature,train_labels, mtry = mtry_best)
      
      train_predicted <- as.data.frame(RF_model_summed$votes)
      train_roc.curve <- roc(train_labels,train_predicted[,which(colnames(train_predicted)==positive)],levels=c(negative,positive),direction=c("<"))
      Summed_signature_auc_in_trainset <- train_roc.curve$auc
    }
  }
  
  #plot signature
  {
    my_comparisons <- list(c("early","late"))
    forplot <- cbind(Fitted_signature_forfit,Summed_signature,Expression_clinical)
    p_fitted <- ggplot(forplot,aes(x=forplot[[clinical_data]],y=Fitted_signature,fill = forplot[[clinical_data]]))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","early"="dark green","late"="black"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Fitted signature: ",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_fitted
    
    p_summed <- ggplot(forplot,aes(x=forplot[[clinical_data]],y=Summed_signature,fill = forplot[[clinical_data]]))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","early"="dark green","late"="black"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Summed signature: ",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_summed
    ggsave(plot = p_fitted,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_training_fitted.pdf"),device = "pdf",height = 5.0,width = 3.8)
    ggsave(plot = p_summed,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_training_summed.pdf"),device = "pdf",height = 5.0,width = 3.8)
  }
  
  
  #validate models by mean AUC and delta AUC on GSE174302 data with early vs. late stage, Gut/exoRBase/GSE174302 data with cancer vs. HD
  {
    TPM_forvalidation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/TPM/GSE174302_TPM_noduplicated_gene.txt",header = TRUE,check.names = FALSE,row.names = 1,sep = "\t")
    clinical_forvalidation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Sample_clinical_information_GSE174302.csv",header = TRUE,check.names = FALSE,row.names = 1)
    rownames(clinical_forvalidation) <- gsub("CRC-2410911","CRC-2410966",rownames(clinical_forvalidation))
    
    TPM_forvalidation <- TPM_forvalidation[,grep("STAD|CRC",colnames(TPM_forvalidation))]
    clinical_forvalidation <- clinical_forvalidation[colnames(TPM_forvalidation),]
    
    
    Marker_genes <- TPM_forvalidation[unique(as.character(Interested_genes$ensembl_gene_id)),]
    rownames(Marker_genes) <- unique(as.character(Interested_genes$ensembl_gene_id))
    Marker_genes[is.na(Marker_genes)] <- 0
    #Marker_genes <- Marker_genes[-grep("NA",rownames(Marker_genes)),]
    Expression_clinical_forvalidation <- cbind(clinical_forvalidation,t(Marker_genes))
    #Expression_clinical_forvalidation <- Expression_clinical_forvalidation[-grep("x",Expression_clinical[[clinical_data]]),]
    
    forvalidation <- Expression_clinical_forvalidation
    
    Summed_signature_forvalidation <- as.data.frame(rowSums(log2(forvalidation[,rownames(Marker_genes)]+1)))
    colnames(Summed_signature_forvalidation) <- "Summed_signature"
    
    #B is a subset of forvalidation, B only contains interested 1 clinical response and multiple genes, colnames same as A
    B <- as.matrix(sapply(forvalidation[,colnames(A)], as.numeric))
    rownames(B) <- rownames(forvalidation)
    Lasso_selected_genes <- as.data.frame(B[,(coefficient@i+1)[-1]])
    genes_coefficients <- coefficient@x[-1]
    
    i=1
    Fitted_signature_forvalidation = 0
    formulas = "Stage_simplified ~ "
    while(i<=ncol(Lasso_selected_genes)){
      Fitted_signature_forvalidation = Fitted_signature_forvalidation + genes_coefficients[i]*Lasso_selected_genes[i]
      formulas <- paste0(formulas,"+",genes_coefficients[i],"*",colnames(Lasso_selected_genes[i]),"")
      i=i+1
    }
    Fitted_signature_forfit <- sigmoid(Fitted_signature_forfit)
    colnames(Fitted_signature_forvalidation) <- "Fitted_signature"
    
    #Fitted score AUC in validation set
    {
      test_labels <- forvalidation[[clinical_data]]
      test_labels <- factor(test_labels, levels = c(negative,positive))
      
      predict_prob <- predict(RF_model_fitted,newdata = Fitted_signature_forvalidation, type = "prob")
    
    #report ROC
    predicted <- predict_prob
    roc.curve <- roc(test_labels,predicted[,which(colnames(predicted)==positive)],levels=c(negative,positive),direction=c("<"))
    ci.auc(roc.curve,conf.level = 0.95)
    record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
    
    Fitted_signature_auc_in_GSE174302 <- roc.curve$auc
    }
    
    #Summed score AUC in validation set
    {
      test_labels <- forvalidation[[clinical_data]]
      test_labels <- factor(test_labels, levels = c(negative,positive))
      
      predict_prob <- predict(RF_model_summed,newdata = Summed_signature_forvalidation, type = "prob")
      
      #report ROC
      predicted <- predict_prob
      roc.curve <- roc(test_labels,predicted[,which(colnames(predicted)==positive)],levels=c(negative,positive),direction=c("<"))
      ci.auc(roc.curve,conf.level = 0.95)
      record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
      
      Summed_signature_auc_in_GSE174302 <- roc.curve$auc
    }
  }
  
  #plot signature
  {
    my_comparisons <- list(c("early","late"))
    forplot <- cbind(Fitted_signature_forvalidation,Summed_signature_forvalidation,clinical_forvalidation)
    p_fitted <- ggplot(forplot,aes(x=forplot[[clinical_data]],y=Fitted_signature,fill = forplot[[clinical_data]]))+
    #p_fitted <- ggplot(forplot,aes(x=forplot$Stage_raw,y=Fitted_signature,fill = forplot[[clinical_data]]))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","early"="dark green","late"="black"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Fitted signature: ",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_fitted
    
    p_summed <- ggplot(forplot,aes(x=forplot[[clinical_data]],y=Summed_signature,fill = forplot[[clinical_data]]))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","early"="dark green","late"="black"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Summed signature: ",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_summed
    ggsave(plot = p_fitted,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_GSE174302_fitted.pdf"),device = "pdf",height = 5.0,width = 3.8)
    ggsave(plot = p_summed,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_GSE174302_summed.pdf"),device = "pdf",height = 5.0,width = 3.8)
  }
  
  message("Training fitted score: ",Fitted_signature_auc_in_trainset)
  message("Training summed score: ",Summed_signature_auc_in_trainset)
  message("Testing Fitted score in GSE174302: ",Fitted_signature_auc_in_GSE174302)
  message("Testing summed score in GSE174302: ",Summed_signature_auc_in_GSE174302)
  result_tmp <- data.frame(gene_signatures[n], formulas, Fitted_signature_auc_in_trainset, Summed_signature_auc_in_trainset, Fitted_signature_auc_in_GSE174302, Summed_signature_auc_in_GSE174302)
  colnames(result_tmp) <- c("Genes","Formula","AUROC of Training fitted score","AUROC of Training summed score","AUROC of Testing Fitted score in GSE174302","AUROC of Testing summed score in GSE174302")
  result <- rbind(result,result_tmp)
  n=n+1
  }
  write.csv(result,paste0(workdir,"/Signature_AUC_Records.csv"))
}

############################################################################################
#Logistic regression to evaluate different gene signatures (continus clinical data)
{
  #variable***
  gene_signature_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Gene_signatures"
  workdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Fitted_signature_forStage_spearman"
  dir.create(workdir)
  gene_signatures <- gsub(".csv","",dir(gene_signature_dir))
  #gene_signatures <- gene_signatures[-grep("Bacterial invasion of epithelial cells_KEGG",gene_signatures)]
  #gene_signatures <- gene_signatures[grep("MolecularCancer",gene_signatures)]
  #gene_signatures <- c("genes_ensembl_gammadeltT_LM22","genes_ensembl_CAF_EPIC")
  
  set.seed(12345)
  n=1
  result <- as.data.frame(matrix(numeric(0),ncol=26))
  colnames(result) <- c("Genes","Formula","Correlation of Training fitted score","Correlation of Training summed score","Correlation of Testing Fitted score in GSE174302","Correlation of Testing summed score in GSE174302","Pvalue of Training fitted score","Pvalue of Training summed score","Pvalue of Testing Fitted score in GSE174302","Pvalue of Testing summed score in GSE174302",
                        "train_wilcox_fitted","train_wilcox_summed","GSE174302_wilcox_fitted","GSE174302_wilcox_summed","GSE133684_wilcox_fitted","GSE133684_wilcox_summed","exoRBase_wilcox_fitted","exoRBase_wilcox_summed",
                        "train_FC_fitted","train_FC_summed","GSE174302_FC_fitted","GSE174302_FC_summed","GSE133684_FC_fitted","GSE133684_FC_summed","exoRBase_FC_fitted","exoRBase_FC_summed")
  while(n<=length(gene_signatures)){
    message(n)
    dir.create(paste0(workdir,"/",gene_signatures[n]))
    setwd(paste0(workdir,"/",gene_signatures[n]))
    
    Interested_genes <- read.csv(paste0(gene_signature_dir,"/",gene_signatures[n],".csv"))
    clinical_data <- "Stage_raw"

    #read in genes of different signatures for training
    {
      TPM_forfit <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/TPM/Plasma_PBMC_tissue_TPM_noduplicated_gene.txt",header = TRUE,check.names = FALSE,row.names = 1,sep = "\t") 
      TPM_forfit <- TPM_forfit[,-grep("NC|PBMC|-T|-N|HCC|LUAD|ESCA|Tumor|Normal|NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",colnames(TPM_forfit))]
      
      clinical <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Sample_clinical_information_x.csv",header = TRUE,check.names = FALSE,row.names = 1)
      rownames(clinical) <- paste0(rownames(clinical),"-pico")
      clinical <- clinical[-grep("NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",rownames(clinical)),]
      clinical <- clinical[colnames(TPM_forfit),]
      
      Marker_genes <- TPM_forfit[unique(as.character(Interested_genes$ensembl_gene_id)),]
      rownames(Marker_genes) <- unique(as.character(Interested_genes$ensembl_gene_id))
      Marker_genes[is.na(Marker_genes)] <- 0 
      #Marker_genes <- Marker_genes[-grep("NA",rownames(Marker_genes)),]
      Expression_clinical <- cbind(clinical,t(Marker_genes))
      Expression_clinical <- Expression_clinical[-grep("x",Expression_clinical[[clinical_data]]),]
      
      forfit <- Expression_clinical
      clinical_forfit <- clinical[rownames(forfit),]
    }
    
    #fit a signature score based on multiomics RNA data, by cv.glmnet
    {
      #A is a subset of forfit, A only contains interested 1 clinical response and multiple genes
      A <- as.matrix(sapply(forfit[,c(clinical_data,unique(as.character(Interested_genes$ensembl_gene_id)))], as.numeric))
      rownames(A) <- rownames(forfit)
      
      cvfit <- cv.glmnet(A[,-1], A[,1] , standardize = TRUE , type.measure = "mse" , nfolds = 5 , alpha = 1) 
      cvfit    
      cvfit$lambda.min
      cvfit$lambda.1se
      coefficient <- coef(cvfit , s = "lambda.min")
      #coefficient <- coef(cvfit , s = "lambda.1se")
      
      if(length((coefficient@i+1)[-1])==0){
        message("Skipped: ",gene_signatures[n],"\nBecause coefficient number = 0")
        n=n+1
        next
      }
      Lasso_selected_genes <- as.data.frame(A[,(coefficient@i+1)[-1]])
      genes_coefficients <- coefficient@x[-1]
      
      i=1
      Fitted_signature_forfit = 0
      formulas = paste0(clinical_data," ~ ")
      while(i<=ncol(Lasso_selected_genes)){
        Fitted_signature_forfit = Fitted_signature_forfit + genes_coefficients[i]*Lasso_selected_genes[i]
        formulas <- paste0(formulas,"+",genes_coefficients[i],"*",colnames(Lasso_selected_genes[i]),"")
        i=i+1
      }
      Fitted_signature_forfit <- sigmoid(Fitted_signature_forfit)
      colnames(Fitted_signature_forfit) <- "Fitted_signature"
    }
    
    #get summed gene set signature
    {
      Summed_signature_forfit <- as.data.frame(rowSums(log2(forfit[,rownames(Marker_genes),drop = FALSE]+1)))
      colnames(Summed_signature_forfit) <- "Summed_signature"
    }
    
    #plot signature
    {
      forplot <- cbind(Fitted_signature_forfit,Summed_signature_forfit,clinical_forfit)
      forplot[[clinical_data]] <- as.numeric(forplot[[clinical_data]])
      p_fitted <-
        ggscatter(forplot, x = "Fitted_signature", y = clinical_data, 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
        xlab("Fitted score")+
        ylab(clinical_data)+
        theme(#legend.position="right",
          legend.position="right",
          panel.grid=element_blank(),
          plot.margin = unit(c(20,20,20,20),"pt"),
          legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.text.x = element_text(face="bold", color="black", size=20, family = "Arial"),
          axis.text.y = element_text(face="bold",  color="black", size=20, family = "Arial"),
          axis.title.x = element_text(face="bold", color="black", size=16, family = "Arial"),
          axis.title.y = element_text(face="bold",color="black", size=24, family = "Arial"))

      p_summed <- 
        ggscatter(forplot, x = "Summed_signature", y = clinical_data, 
                            add = "reg.line", conf.int = TRUE, 
                            cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
        xlab("Summed score")+
        ylab(clinical_data)+
        theme(#legend.position="right",
          legend.position="right",
          panel.grid=element_blank(),
          plot.margin = unit(c(20,20,20,20),"pt"),
          legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.text.x = element_text(face="bold", color="black", size=20, family = "Arial"),
          axis.text.y = element_text(face="bold",  color="black", size=20, family = "Arial"),
          axis.title.x = element_text(face="bold", color="black", size=16, family = "Arial"),
          axis.title.y = element_text(face="bold",color="black", size=24, family = "Arial"))

      ggsave(plot = p_fitted,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_training_fitted_spearman.pdf"),device = "pdf",height = 5.0,width = 5)
      ggsave(plot = p_summed,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_training_summed_spearman.pdf"),device = "pdf",height = 5.0,width = 5)
    }
    
    
    #train summed signature and fitted score seperately by spearman correlation
    {
      #Fitted score R in train set by spearman correlation
      {
        R_fitted <- spearman.test(forplot$Fitted_signature,as.numeric(forplot[[clinical_data]]),alternative = "two.sided", approximation = "exact")
        Fitted_signature_R_in_trainset <- R_fitted$estimate
        Fitted_signature_Pvalue_in_trainset <- R_fitted$p.value
      }
      
      #Summed score R in train set by spearman correlation
      {
        R_Summed <- spearman.test(forplot$Summed_signature,as.numeric(forplot[[clinical_data]]),alternative = "two.sided", approximation = "exact")
        Summed_signature_R_in_trainset <- R_Summed$estimate
        Summed_signature_Pvalue_in_trainset <- R_Summed$p.value
      }
    }
    
    #validate models by correlation on GSE174302 data
    {
      TPM_forvalidation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/TPM/GSE174302_TPM_noduplicated_gene.txt",header = TRUE,check.names = FALSE,row.names = 1,sep = "\t")
      clinical_forvalidation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Sample_clinical_information_GSE174302.csv",header = TRUE,check.names = FALSE,row.names = 1)
      rownames(clinical_forvalidation) <- gsub("CRC-2410911","CRC-2410966",rownames(clinical_forvalidation))
      
      TPM_forvalidation <- TPM_forvalidation[,grep("STAD|CRC",colnames(TPM_forvalidation))]
      clinical_forvalidation <- clinical_forvalidation[colnames(TPM_forvalidation),]
      
      
      Marker_genes <- TPM_forvalidation[unique(as.character(Interested_genes$ensembl_gene_id)),]
      rownames(Marker_genes) <- unique(as.character(Interested_genes$ensembl_gene_id))
      Marker_genes[is.na(Marker_genes)] <- 0
      #Marker_genes <- Marker_genes[-grep("NA",rownames(Marker_genes)),]
      Expression_clinical_forvalidation <- cbind(clinical_forvalidation,t(Marker_genes))
      #Expression_clinical_forvalidation <- Expression_clinical_forvalidation[-grep("x",Expression_clinical[[clinical_data]]),]
      
      forvalidation <- Expression_clinical_forvalidation
      
      Summed_signature_forvalidation <- as.data.frame(rowSums(log2(forvalidation[,rownames(Marker_genes),drop = FALSE]+1)))
      colnames(Summed_signature_forvalidation) <- "Summed_signature"
      
      #B is a subset of forvalidation, B only contains interested 1 clinical response and multiple genes, colnames same as A
      B <- as.matrix(sapply(forvalidation[,colnames(A)], as.numeric))
      rownames(B) <- rownames(forvalidation)
      Lasso_selected_genes <- as.data.frame(B[,(coefficient@i+1)[-1]])
      genes_coefficients <- coefficient@x[-1]
      
      i=1
      Fitted_signature_forvalidation = 0
      formulas = paste0(clinical_data," ~ ")
      while(i<=ncol(Lasso_selected_genes)){
        Fitted_signature_forvalidation = Fitted_signature_forvalidation + genes_coefficients[i]*Lasso_selected_genes[i]
        formulas <- paste0(formulas,"+",genes_coefficients[i],"*",colnames(Lasso_selected_genes[i]),"")
        i=i+1
      }
      Fitted_signature_forfit <- sigmoid(Fitted_signature_forfit)
      colnames(Fitted_signature_forvalidation) <- "Fitted_signature"
    }
    
    #plot signature
    {
      forplot <- cbind(Fitted_signature_forvalidation,Summed_signature_forvalidation,clinical_forvalidation)
      forplot[[clinical_data]] <- as.numeric(forplot[[clinical_data]])
      p_fitted <-
        ggscatter(forplot, x = "Fitted_signature", y = clinical_data, 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
        xlab("Fitted score")+
        ylab(clinical_data)+
        theme(#legend.position="right",
          legend.position="right",
          panel.grid=element_blank(),
          plot.margin = unit(c(20,20,20,20),"pt"),
          legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.text.x = element_text(face="bold", color="black", size=20, family = "Arial"),
          axis.text.y = element_text(face="bold",  color="black", size=20, family = "Arial"),
          axis.title.x = element_text(face="bold", color="black", size=16, family = "Arial"),
          axis.title.y = element_text(face="bold",color="black", size=24, family = "Arial"))
      
      p_summed <- 
        ggscatter(forplot, x = "Summed_signature", y = clinical_data, 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
        xlab("Summed score")+
        ylab(clinical_data)+
        theme(#legend.position="right",
          legend.position="right",
          panel.grid=element_blank(),
          plot.margin = unit(c(20,20,20,20),"pt"),
          legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.text.x = element_text(face="bold", color="black", size=20, family = "Arial"),
          axis.text.y = element_text(face="bold",  color="black", size=20, family = "Arial"),
          axis.title.x = element_text(face="bold", color="black", size=16, family = "Arial"),
          axis.title.y = element_text(face="bold",color="black", size=24, family = "Arial"))

      ggsave(plot = p_fitted,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_GSE174302_fitted_spearman.pdf"),device = "pdf",height = 5.0,width = 3.8)
      ggsave(plot = p_summed,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_GSE174302_summed_spearman.pdf"),device = "pdf",height = 5.0,width = 3.8)
    }
    
    #validation summed signature and fitted score seperately by spearman correlation
    {
    #Fitted score AUC in validation set
    {
      R_GSE174302 <- spearman.test(forplot$Fitted_signature,as.numeric(forplot[[clinical_data]]),alternative = "two.sided", approximation = "exact")
      Fitted_signature_R_in_GSE174302 <- R_GSE174302$estimate
      Fitted_signature_Pvalue_in_GSE174302 <- R_GSE174302$p.value
    }
    
    #Summed score AUC in validation set
    {
      R_GSE174302 <- spearman.test(forplot$Summed_signature,as.numeric(forplot[[clinical_data]]),alternative = "two.sided", approximation = "exact")
      Summed_signature_R_in_GSE174302 <- R_GSE174302$estimate
      Summed_signature_Pvalue_in_GSE174302 <- R_GSE174302$p.value
    }
    }
    
    message("Training fitted score: ",Fitted_signature_R_in_trainset)
    message("Training summed score: ",Summed_signature_R_in_trainset)
    message("Testing Fitted score in GSE174302: ",Fitted_signature_R_in_GSE174302)
    message("Testing summed score in GSE174302: ",Summed_signature_R_in_GSE174302)
    result_tmp <- data.frame(gene_signatures[n], formulas, Fitted_signature_R_in_trainset, Summed_signature_R_in_trainset, Fitted_signature_R_in_GSE174302, Summed_signature_R_in_GSE174302,Fitted_signature_Pvalue_in_trainset, Summed_signature_Pvalue_in_trainset, Fitted_signature_Pvalue_in_GSE174302, Summed_signature_Pvalue_in_GSE174302)
    colnames(result_tmp) <- c("Genes","Formula","Correlation of Training fitted score","Correlation of Training summed score","Correlation of Testing Fitted score in GSE174302","Correlation of Testing summed score in GSE174302","Pvalue of Training fitted score","Pvalue of Training summed score","Pvalue of Testing Fitted score in GSE174302","Pvalue of Testing summed score in GSE174302")
    result <- rbind(result,result_tmp)
    n=n+1
  }
  write.csv(result,paste0(workdir,"/Signature_AUC_Records.csv"))
}

#########################################################################
#Logistic regression to evaluate different gene signatures in GI vs. HD
{
  #KEGG_pathways <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/09.Published_data/pathway/PATH_ID_NAME_modified.csv",header = TRUE)
  #KEGG_pathways$DESCRPTION <- gsub("/","-",KEGG_pathways$DESCRPTION)
  #pathway_description <- unique(KEGG_pathways$DESCRPTION)
  #i=1
  #while(i<=length(pathway_description)){
  #  single_pathway <- data.frame("Gene_Name"=KEGG_pathways[which(KEGG_pathways$DESCRPTION==pathway_description[i]),]$hgnc_symbol,"ensembl_gene_id"=KEGG_pathways[which(KEGG_pathways$DESCRPTION==pathway_description[i]),]$ensembl_gene_id)
  #  write.csv(single_pathway,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Gene_signatures/",pathway_description[i],"_KEGG.csv"),quote = FALSE)
  #  i=i+1
  #}  
  
  
  #variable***
  gene_signature_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Gene_signatures"
  workdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Fitted_signature_forDiagnosis"
  dir.create(workdir)
  gene_signatures <- gsub(".csv","",dir(gene_signature_dir))
  #gene_signatures <- c("B cell receptor signaling pathway")
  
  n=1
  result <- as.data.frame(matrix(numeric(0),ncol=10))
  colnames(result) <- c("Genes","Formula","AUROC of training fitted score","AUROC of training summed score","AUROC of testing fitted score in GSE174302","AUROC of testing summed score in GSE174302","AUROC of testing fitted score in GSE133684","AUROC of testing summed score in GSE133684","AUROC of testing fitted score in exoRBase","AUROC of testing summed score in exoRBase")
  while(n<=length(gene_signatures)){
    
    dir.create(paste0(workdir,"/",gene_signatures[n]))
    setwd(paste0(workdir,"/",gene_signatures[n]))
    
  Interested_genes <- read.csv(paste0(gene_signature_dir,"/",gene_signatures[n],".csv"))
  #Interested_genes <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/B cell receptor signaling pathway.csv")
  clinical_data <- "Label" 
  negative <- "Healthy" 
  positive <- "Cancer"
  
  #read in genes of different signatures for training
  {
    TPM_forfit <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/TPM/Plasma_PBMC_tissue_TPM_noduplicated_gene.txt",header = TRUE,check.names = FALSE,row.names = 1,sep = "\t") 
    TPM_forfit <- TPM_forfit[,-grep("PBMC|-T|-N|Tumor|Normal|NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",colnames(TPM_forfit))]
    
    clinical <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Multiomics_sample_label.csv",header = TRUE,check.names = FALSE,row.names = 1)
    #clinical <- clinical[-grep("NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",rownames(clinical)),]
    clinical <- clinical[colnames(TPM_forfit),,drop = FALSE]
    
    Marker_genes <- TPM_forfit[unique(as.character(Interested_genes$ensembl_gene_id)),]
    rownames(Marker_genes) <- unique(as.character(Interested_genes$ensembl_gene_id))
    Marker_genes[is.na(Marker_genes)] <- 0
    #Marker_genes <- Marker_genes[-grep("NA",rownames(Marker_genes)),]
    Expression_clinical <- cbind(clinical,t(Marker_genes))
  
    forfit <- Expression_clinical
  }
  
  #fit a signature score based on multiomics RNA data, by cv.glmnet
  {
    #A is a subset of forfit, A only contains interested 1 clinical response and multiple genes
    A <- as.matrix(sapply(forfit[,c(clinical_data,unique(as.character(Interested_genes$ensembl_gene_id)))], as.numeric))
    rownames(A) <- rownames(forfit)
    
    cvfit <- cv.glmnet(A[,-1], A[,1] , standardize = TRUE , type.measure = "auc", family = "binomial" , nfolds = 5, alpha = 1) 
    cvfit    
    cvfit$lambda.min
    cvfit$lambda.1se
    coefficient <- coef(cvfit , s = "lambda.min")
    #coefficient <- coef(cvfit , s = "lambda.1se")
    
    Lasso_selected_genes <- as.data.frame(A[,(coefficient@i+1)[-1]])
    genes_coefficients <- coefficient@x[-1]
    
    i=1
    Fitted_signature_forfit = 0
    formulas = paste0(clinical_data," ~ ")
    while(i<=ncol(Lasso_selected_genes)){
      Fitted_signature_forfit = Fitted_signature_forfit + genes_coefficients[i]*Lasso_selected_genes[i]
      formulas <- paste0(formulas,"+",genes_coefficients[i],"*",colnames(Lasso_selected_genes[i]),"")
      i=i+1
    }
    Fitted_signature_forfit <- sigmoid(Fitted_signature_forfit)
    colnames(Fitted_signature_forfit) <- "Fitted_signature"
  }
  
  #train summed signature and fitted score seperately by caret
  {
    library(caret)
    library(pROC)
    #Fitted score AUC in train set by tree model
    {
      train_labels <- forfit[[clinical_data]]
      train_labels <- factor(train_labels, levels = c(negative,positive))
      set.seed(100000)
      #fix model parameters, such as mtry, ntree
      tune_mtry <- tuneRF(Fitted_signature_forfit,train_labels) #tune mtry
      tune_mtry <- as.data.frame(tune_mtry)
      mtry_best <- tune_mtry[order(tune_mtry$OOBError,decreasing=FALSE),][1,1]
      
      #train RF model
      RF_model_fitted <- randomForest(Fitted_signature_forfit,train_labels, mtry = mtry_best)
      
      train_predicted <- as.data.frame(RF_model_fitted$votes)
      train_roc.curve <- roc(train_labels,train_predicted[,which(colnames(train_predicted)==positive)],levels=c(negative,positive),direction=c("<"))
      Fitted_signature_auc_in_trainset <- train_roc.curve$auc
    }
    
    #Fitted score AUC in train set by LR model
    #There must be two or more independent variables, or predictors, for a logistic regression.
    {
      #params.grid <- expand.grid(alpha = c(0,0.5,1),lambda = c(0,0.01,0.1,1))
      # alpha: relative weighting of L1 and L2 regularization
      # lambda: degree of regularization, see glmnet documentation for detail
      
      # train是caret封装的一个用于调参的函数
      # 和rfeControl类似，trainControl也返回一个列表，用于控制train函数的行为
      # method同样可以是"boot","cv"等，number参数控制了交叉验证的折数和bootstraping的重复次数
      #tr.ctrl <- trainControl(method="cv",number = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
      #cv.fitted <- train(Fitted_signature_forfit,train_labels,method="glmnet",family="binomial",metric = "ROC",tuneGrid = params.grid,preProcess = NULL,trControl = tr.ctrl )
    }
    
    
    #get summed gene set signature
    {
      Summed_signature <- as.data.frame(rowSums(log2(forfit[,rownames(Marker_genes)]+1)))
      colnames(Summed_signature) <- "Summed_signature"
    }
    #Summed score AUC in train set
    {
      train_labels <- forfit[[clinical_data]]
      train_labels <- factor(train_labels, levels = c(negative,positive))
      set.seed(100000)
      #fix model parameters, such as mtry, ntree
      tune_mtry <- tuneRF(Summed_signature,train_labels) #tune mtry
      tune_mtry <- as.data.frame(tune_mtry)
      mtry_best <- tune_mtry[order(tune_mtry$OOBError,decreasing=FALSE),][1,1]
      
      #train RF model use n_sample-1
      RF_model_summed <- randomForest(Summed_signature,train_labels, mtry = mtry_best)
      
      train_predicted <- as.data.frame(RF_model_summed$votes)
      train_roc.curve <- roc(train_labels,train_predicted[,which(colnames(train_predicted)==positive)],levels=c(negative,positive),direction=c("<"))
      Summed_signature_auc_in_trainset <- train_roc.curve$auc
    }
  }
  #plot signature
  {
    my_comparisons <- list(c("Cancer","Healthy"))
    forplot <- cbind(Fitted_signature_forfit,Summed_signature,clinical)
    p_fitted <- ggplot(forplot,aes(x=Label,y=Fitted_signature,fill = Label))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","Cancer"="red","Healthy"="blue"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Fitted signature: ",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_fitted
    
    p_summed <- ggplot(forplot,aes(x=Label,y=Summed_signature,fill = Label))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","Cancer"="red","Healthy"="blue"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Summed signature: ",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_summed
    
    wilcox <- wilcox.test(forplot[forplot$Label=="Cancer",]$Fitted_signature,forplot[forplot$Label=="Healthy",]$Fitted_signature)
    train_wilcox_fitted <- wilcox$p.value
    train_FC_fitted <- mean(forplot[forplot$Label=="Cancer",]$Fitted_signature)/mean(forplot[forplot$Label=="Healthy",]$Fitted_signature)
    wilcox <- wilcox.test(forplot[forplot$Label=="Cancer",]$Summed_signature,forplot[forplot$Label=="Healthy",]$Summed_signature)
    train_wilcox_summed <- wilcox$p.value
    train_FC_summed <- mean(forplot[forplot$Label=="Cancer",]$Summed_signature)/mean(forplot[forplot$Label=="Healthy",]$Summed_signature)
    ggsave(plot = p_fitted,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_training_fitted.pdf"),device = "pdf",height = 5.0,width = 3.8)
    ggsave(plot = p_summed,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_training_summed.pdf"),device = "pdf",height = 5.0,width = 3.8)
  }
  
  
  #validate models by mean AUC and delta AUC on GSE174302 data with early vs. late stage, Gut/exoRBase/GSE174302 data with cancer vs. HD, and generate boxplot
  #GSE174302
  #AUC
  {
    TPM_forvalidation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/TPM/GSE174302_TPM_noduplicated_gene.txt",header = TRUE,check.names = FALSE,row.names = 1,sep = "\t")
    clinical_forvalidation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/GSE174302_sample_label.csv",header = TRUE,check.names = FALSE,row.names = 1)
    
    #TPM_forvalidation <- TPM_forvalidation[,grep("STAD|CRC|NC",colnames(TPM_forvalidation))]
    clinical_forvalidation <- clinical_forvalidation[colnames(TPM_forvalidation),,drop = FALSE]
    
    Marker_genes <- TPM_forvalidation[unique(as.character(Interested_genes$ensembl_gene_id)),]
    rownames(Marker_genes) <- unique(as.character(Interested_genes$ensembl_gene_id))
    Marker_genes[is.na(Marker_genes)] <- 0
    #Marker_genes <- Marker_genes[-grep("NA",rownames(Marker_genes)),]
    Expression_clinical_forvalidation <- cbind(clinical_forvalidation,t(Marker_genes))
    #Expression_clinical_forvalidation <- Expression_clinical_forvalidation[-which(Expression_clinical_forvalidation$Stage_simplified=="x"),]
    
    forvalidation <- Expression_clinical_forvalidation
    
    Summed_signature_forvalidation <- as.data.frame(rowSums(log2(forvalidation[,rownames(Marker_genes)]+1)))
    colnames(Summed_signature_forvalidation) <- "Summed_signature"
    
    #B is a subset of forvalidation, B only contains interested 1 clinical response and multiple genes, colnames same as A
    B <- as.matrix(sapply(forvalidation[,colnames(A)], as.numeric))
    rownames(B) <- rownames(forvalidation)
    Lasso_selected_genes <- as.data.frame(B[,(coefficient@i+1)[-1]])
    genes_coefficients <- coefficient@x[-1]
    
    i=1
    Fitted_signature_forvalidation = 0
    formulas = paste0(clinical_data, " ~ ")
    while(i<=ncol(Lasso_selected_genes)){
      Fitted_signature_forvalidation = Fitted_signature_forvalidation + genes_coefficients[i]*Lasso_selected_genes[i]
      formulas <- paste0(formulas,"+",genes_coefficients[i],"*",colnames(Lasso_selected_genes[i]),"")
      i=i+1
    }
    Fitted_signature_forvalidation <- sigmoid(Fitted_signature_forvalidation)
    colnames(Fitted_signature_forvalidation) <- "Fitted_signature"
    
    #Fitted score AUC in validation set
    {
      test_labels <- forvalidation[[clinical_data]]
      test_labels <- factor(test_labels, levels = c(negative,positive))
      
      predict_prob <- predict(RF_model_fitted,newdata = Fitted_signature_forvalidation, type = "prob")
      
      #report ROC
      predicted <- predict_prob
      roc.curve <- roc(test_labels,predicted[,which(colnames(predicted)==positive)],levels=c(negative,positive),direction=c("<"))
      ci.auc(roc.curve,conf.level = 0.95)
      record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
      
      Fitted_signature_auc_in_GSE174302 <- roc.curve$auc
    }
    
    #Summed score AUC in validation set
    {
      test_labels <- forvalidation[[clinical_data]]
      test_labels <- factor(test_labels, levels = c(negative,positive))
      
      predict_prob <- predict(RF_model_summed,newdata = Summed_signature_forvalidation, type = "prob")
      
      #report ROC
      predicted <- predict_prob
      roc.curve <- roc(test_labels,predicted[,which(colnames(predicted)==positive)],levels=c(negative,positive),direction=c("<"))
      ci.auc(roc.curve,conf.level = 0.95)
      record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
      
      Summed_signature_auc_in_GSE174302 <- roc.curve$auc
    }
  }
  #plot signature
  {
    my_comparisons <- list(c("Cancer","Healthy"))
    forplot <- cbind(Fitted_signature_forvalidation,Summed_signature_forvalidation,clinical_forvalidation)
    p_fitted <- ggplot(forplot,aes(x=Label,y=Fitted_signature,fill = Label))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","Cancer"="red","Healthy"="blue"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Fitted signature",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_fitted
    
    p_summed <- ggplot(forplot,aes(x=Label,y=Summed_signature,fill = Label))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","Cancer"="red","Healthy"="blue"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Summed signature: ",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_summed
    wilcox <- wilcox.test(forplot[forplot$Label=="Cancer",]$Fitted_signature,forplot[forplot$Label=="Healthy",]$Fitted_signature)
    GSE174302_wilcox_fitted <- wilcox$p.value
    GSE174302_FC_fitted <- mean(forplot[forplot$Label=="Cancer",]$Fitted_signature)/mean(forplot[forplot$Label=="Healthy",]$Fitted_signature)
    wilcox <- wilcox.test(forplot[forplot$Label=="Cancer",]$Summed_signature,forplot[forplot$Label=="Healthy",]$Summed_signature)
    GSE174302_wilcox_summed <- wilcox$p.value
    GSE174302_FC_summed <- mean(forplot[forplot$Label=="Cancer",]$Summed_signature)/mean(forplot[forplot$Label=="Healthy",]$Summed_signature)
    ggsave(plot = p_fitted,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_testing_GSE174302_fitted.pdf"),device = "pdf",height = 5.0,width = 3.8)
    ggsave(plot = p_summed,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_testing_GSE174302_summed.pdf"),device = "pdf",height = 5.0,width = 3.8)
  }
  #Gut
  {
    TPM_forvalidation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/GSE133684_count_matrix_TPM.txt",header = TRUE,check.names = FALSE,row.names = 1,sep = "\t")
    clinical_forvalidation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/GSE133684_sample_label.csv",header = TRUE,check.names = FALSE,row.names = 1)
    
    #TPM_forvalidation <- TPM_forvalidation[,grep("STAD|CRC|NC",colnames(TPM_forvalidation))]
    TPM_forvalidation <- TPM_forvalidation[,rownames(clinical_forvalidation),drop = FALSE]
    
    Marker_genes <- TPM_forvalidation[unique(as.character(Interested_genes$ensembl_gene_id)),]
    rownames(Marker_genes) <- unique(as.character(Interested_genes$ensembl_gene_id))
    Marker_genes[is.na(Marker_genes)] <- 0
    #Marker_genes <- Marker_genes[-grep("NA",rownames(Marker_genes)),]
    Expression_clinical_forvalidation <- cbind(clinical_forvalidation,t(Marker_genes))
    #Expression_clinical_forvalidation <- Expression_clinical_forvalidation[-which(Expression_clinical_forvalidation$Stage_simplified=="x"),]
    
    forvalidation <- Expression_clinical_forvalidation
    
    Summed_signature_forvalidation <- as.data.frame(rowSums(log2(forvalidation[,rownames(Marker_genes)]+1)))
    colnames(Summed_signature_forvalidation) <- "Summed_signature"
    
    #B is a subset of forvalidation, B only contains interested 1 clinical response and multiple genes, colnames same as A
    B <- as.matrix(sapply(forvalidation[,colnames(A)], as.numeric))
    rownames(B) <- rownames(forvalidation)
    Lasso_selected_genes <- as.data.frame(B[,(coefficient@i+1)[-1]])
    genes_coefficients <- coefficient@x[-1]
    
    i=1
    Fitted_signature_forvalidation = 0
    formulas = paste0(clinical_data, " ~ ")
    while(i<=ncol(Lasso_selected_genes)){
      Fitted_signature_forvalidation = Fitted_signature_forvalidation + genes_coefficients[i]*Lasso_selected_genes[i]
      formulas <- paste0(formulas,"+",genes_coefficients[i],"*",colnames(Lasso_selected_genes[i]),"")
      i=i+1
    }
    Fitted_signature_forvalidation <- sigmoid(Fitted_signature_forvalidation)
    colnames(Fitted_signature_forvalidation) <- "Fitted_signature"
    
    #Fitted score AUC in validation set
    {
      test_labels <- forvalidation[[clinical_data]]
      test_labels <- factor(test_labels, levels = c(negative,positive))
      
      predict_prob <- predict(RF_model_fitted,newdata = Fitted_signature_forvalidation, type = "prob")
      
      #report ROC
      predicted <- predict_prob
      roc.curve <- roc(test_labels,predicted[,which(colnames(predicted)==positive)],levels=c(negative,positive),direction=c("<"))
      ci.auc(roc.curve,conf.level = 0.95)
      record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
      
      Fitted_signature_auc_in_GSE133684 <- roc.curve$auc
    }
    
    #Summed score AUC in validation set
    {
      test_labels <- forvalidation[[clinical_data]]
      test_labels <- factor(test_labels, levels = c(negative,positive))
      
      predict_prob <- predict(RF_model_summed,newdata = Summed_signature_forvalidation, type = "prob")
      
      #report ROC
      predicted <- predict_prob
      roc.curve <- roc(test_labels,predicted[,which(colnames(predicted)==positive)],levels=c(negative,positive),direction=c("<"))
      ci.auc(roc.curve,conf.level = 0.95)
      record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
      
      Summed_signature_auc_in_GSE133684 <- roc.curve$auc
    }
  }
  #plot signature
  {
    my_comparisons <- list(c("Cancer","Healthy"))
    forplot <- cbind(Fitted_signature_forvalidation,Summed_signature_forvalidation,clinical_forvalidation)
    p_fitted <- ggplot(forplot,aes(x=Label,y=Fitted_signature,fill = Label))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","Cancer"="red","Healthy"="blue"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Fitted signature: ",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_fitted
    
    p_summed <- ggplot(forplot,aes(x=Label,y=Summed_signature,fill = Label))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","Cancer"="red","Healthy"="blue"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Summed signature: ",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_summed
    wilcox <- wilcox.test(forplot[forplot$Label=="Cancer",]$Fitted_signature,forplot[forplot$Label=="Healthy",]$Fitted_signature)
    GSE133684_wilcox_fitted <- wilcox$p.value
    GSE133684_FC_fitted <- mean(forplot[forplot$Label=="Cancer",]$Fitted_signature)/mean(forplot[forplot$Label=="Healthy",]$Fitted_signature)
    wilcox <- wilcox.test(forplot[forplot$Label=="Cancer",]$Summed_signature,forplot[forplot$Label=="Healthy",]$Summed_signature)
    GSE133684_wilcox_summed <- wilcox$p.value
    GSE133684_FC_summed <- mean(forplot[forplot$Label=="Cancer",]$Summed_signature)/mean(forplot[forplot$Label=="Healthy",]$Summed_signature)
    ggsave(plot = p_fitted,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_testing_GSE133684_fitted.pdf"),device = "pdf",height = 5.0,width = 3.8)
    ggsave(plot = p_summed,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_testing_GSE133684_summed.pdf"),device = "pdf",height = 5.0,width = 3.8)
  }
  #exoRBase
  {
    TPM_forvalidation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/TPM/exoRbase_featurecounts_inhouse_TPM.txt",header = TRUE,check.names = FALSE,row.names = 1,sep = "\t")
    clinical_forvalidation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/exoRBase_sample_label.csv",header = TRUE,check.names = FALSE,row.names = 1)
    
    #TPM_forvalidation <- TPM_forvalidation[,grep("STAD|CRC|NC",colnames(TPM_forvalidation))]
    TPM_forvalidation <- TPM_forvalidation[,rownames(clinical_forvalidation),drop = FALSE]
    
    Marker_genes <- TPM_forvalidation[unique(as.character(Interested_genes$ensembl_gene_id)),]
    rownames(Marker_genes) <- unique(as.character(Interested_genes$ensembl_gene_id))
    Marker_genes[is.na(Marker_genes)] <- 0
    #Marker_genes <- Marker_genes[-grep("NA",rownames(Marker_genes)),]
    Expression_clinical_forvalidation <- cbind(clinical_forvalidation,t(Marker_genes))
    #Expression_clinical_forvalidation <- Expression_clinical_forvalidation[-which(Expression_clinical_forvalidation$Stage_simplified=="x"),]
    
    forvalidation <- Expression_clinical_forvalidation
    
    Summed_signature_forvalidation <- as.data.frame(rowSums(log2(forvalidation[,rownames(Marker_genes)]+1)))
    colnames(Summed_signature_forvalidation) <- "Summed_signature"
    
    #B is a subset of forvalidation, B only contains interested 1 clinical response and multiple genes, colnames same as A
    B <- as.matrix(sapply(forvalidation[,colnames(A)], as.numeric))
    rownames(B) <- rownames(forvalidation)
    Lasso_selected_genes <- as.data.frame(B[,(coefficient@i+1)[-1]])
    genes_coefficients <- coefficient@x[-1]
    
    i=1
    Fitted_signature_forvalidation = 0
    formulas = paste0(clinical_data, " ~ ")
    while(i<=ncol(Lasso_selected_genes)){
      Fitted_signature_forvalidation = Fitted_signature_forvalidation + genes_coefficients[i]*Lasso_selected_genes[i]
      formulas <- paste0(formulas,"+",genes_coefficients[i],"*",colnames(Lasso_selected_genes[i]),"")
      i=i+1
    }
    Fitted_signature_forvalidation <- sigmoid(Fitted_signature_forvalidation)
    colnames(Fitted_signature_forvalidation) <- "Fitted_signature"
    
    #Fitted score AUC in validation set
    {
      test_labels <- forvalidation[[clinical_data]]
      test_labels <- factor(test_labels, levels = c(negative,positive))
      
      predict_prob <- predict(RF_model_fitted,newdata = Fitted_signature_forvalidation, type = "prob")
      
      #report ROC
      predicted <- predict_prob
      roc.curve <- roc(test_labels,predicted[,which(colnames(predicted)==positive)],levels=c(negative,positive),direction=c("<"))
      ci.auc(roc.curve,conf.level = 0.95)
      record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
      
      Fitted_signature_auc_in_exoRBase <- roc.curve$auc
    }
    
    #Summed score AUC in validation set
    {
      test_labels <- forvalidation[[clinical_data]]
      test_labels <- factor(test_labels, levels = c(negative,positive))
      
      predict_prob <- predict(RF_model_summed,newdata = Summed_signature_forvalidation, type = "prob")
      
      #report ROC
      predicted <- predict_prob
      roc.curve <- roc(test_labels,predicted[,which(colnames(predicted)==positive)],levels=c(negative,positive),direction=c("<"))
      ci.auc(roc.curve,conf.level = 0.95)
      record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
      
      Summed_signature_auc_in_exoRBase <- roc.curve$auc
    }
  }
  #plot signature
  {
    my_comparisons <- list(c("Cancer","Healthy"))
    forplot <- cbind(Fitted_signature_forvalidation,Summed_signature_forvalidation,clinical_forvalidation)
    p_fitted <- ggplot(forplot,aes(x=Label,y=Fitted_signature,fill = Label))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","Cancer"="red","Healthy"="blue"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Fitted signature: ",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_fitted
    
    p_summed <- ggplot(forplot,aes(x=Label,y=Summed_signature,fill = Label))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
      scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","Cancer"="red","Healthy"="blue"))+
      #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      xlab("")+
      ylab(paste0("Summed signature: ",gene_signatures[n]))+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(face="bold",  color="black", size=24),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         method.args = list(alternative = "two.sided",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p_summed
    wilcox <- wilcox.test(forplot[forplot$Label=="Cancer",]$Fitted_signature,forplot[forplot$Label=="Healthy",]$Fitted_signature)
    exoRBase_wilcox_fitted <- wilcox$p.value
    exoRBase_FC_fitted <- mean(forplot[forplot$Label=="Cancer",]$Fitted_signature)/mean(forplot[forplot$Label=="Healthy",]$Fitted_signature)
    wilcox <- wilcox.test(forplot[forplot$Label=="Cancer",]$Summed_signature,forplot[forplot$Label=="Healthy",]$Summed_signature)
    exoRBase_wilcox_summed <- wilcox$p.value
    exoRBase_FC_summed <- mean(forplot[forplot$Label=="Cancer",]$Summed_signature)/mean(forplot[forplot$Label=="Healthy",]$Summed_signature)
    ggsave(plot = p_fitted,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_testing_exoRBase_fitted.pdf"),device = "pdf",height = 5.0,width = 3.8)
    ggsave(plot = p_summed,path = paste0(workdir,"/",gene_signatures[n],"/"), filename = paste0(gene_signatures[n],"_testing_exoRBase_summed.pdf"),device = "pdf",height = 5.0,width = 3.8)
  }
  
  message("Training fitted score: ",Fitted_signature_auc_in_trainset)
  message("Training summed score: ",Summed_signature_auc_in_trainset)
  message("Testing Fitted score in GSE174302: ",Fitted_signature_auc_in_GSE174302)
  message("Testing summed score in GSE174302: ",Summed_signature_auc_in_GSE174302)
  message("Testing Fitted score in GSE133684: ",Fitted_signature_auc_in_GSE133684)
  message("Testing summed score in GSE133684: ",Summed_signature_auc_in_GSE133684)
  message("Testing Fitted score in exoRBase: ",Fitted_signature_auc_in_exoRBase)
  message("Testing summed score in exoRBase: ",Summed_signature_auc_in_exoRBase)
  
  result_tmp <- data.frame(gene_signatures[n], formulas, Fitted_signature_auc_in_trainset, Summed_signature_auc_in_trainset, Fitted_signature_auc_in_GSE174302, Summed_signature_auc_in_GSE174302,Fitted_signature_auc_in_GSE133684,Summed_signature_auc_in_GSE133684,Fitted_signature_auc_in_exoRBase,Summed_signature_auc_in_exoRBase,
                           train_wilcox_fitted,train_wilcox_summed,GSE174302_wilcox_fitted,GSE174302_wilcox_summed,GSE133684_wilcox_fitted,GSE133684_wilcox_summed,exoRBase_wilcox_fitted,exoRBase_wilcox_summed,
                           train_FC_fitted,train_FC_summed,GSE174302_FC_fitted,GSE174302_FC_summed,GSE133684_FC_fitted,GSE133684_FC_summed,exoRBase_FC_fitted,exoRBase_FC_summed)
  colnames(result_tmp) <- c("Genes","Formula","AUROC of training fitted score","AUROC of training summed score","AUROC of testing fitted score in GSE174302","AUROC of testing summed score in GSE174302","AUROC of testing fitted score in GSE133684","AUROC of testing summed score in GSE133684","AUROC of testing fitted score in exoRBase","AUROC of testing summed score in exoRBase",
                            "train_wilcox_fitted","train_wilcox_summed","GSE174302_wilcox_fitted","GSE174302_wilcox_summed","GSE133684_wilcox_fitted","GSE133684_wilcox_summed","exoRBase_wilcox_fitted","exoRBase_wilcox_summed",
                            "train_FC_fitted","train_FC_summed","GSE174302_FC_fitted","GSE174302_FC_summed","GSE133684_FC_fitted","GSE133684_FC_summed","exoRBase_FC_fitted","exoRBase_FC_summed")
  result <- rbind(result,result_tmp)
  n=n+1
  }
  write.csv(result,paste0(workdir,"/Signature_AUC_Records.csv"))
}

#########################################################################
#Paired correlation between tissue and plasma (gene level)
{
  Plasma <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Plasma_TPM.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  Tissue <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Tissue_TPM.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  
  gene_ids <- rownames(Plasma)
  
  sample_ids <- c("CRC-PKU-27","CRC-PKU-28","CRC-PKU-29","CRC-PKU-30","CRC-PKU-32","CRC-PKU-34",
                  "CRC-PKU-35","CRC-PKU-36","CRC-PKU-37","CRC-PKU-38","CRC-PKU-39","CRC-PKU-40","CRC-PKU-41")
  
  Tissue_T_forcor <- Tissue[gene_ids,paste0(sample_ids,"-T")]
  colnames(Tissue_T_forcor) <- gsub("-T","",colnames(Tissue_T_forcor))
  Tissue_N_forcor <- Tissue[gene_ids,paste0(sample_ids,"-N")]
  colnames(Tissue_N_forcor) <- gsub("-N","",colnames(Tissue_N_forcor))
  Plasma_forcor <- Plasma[gene_ids,paste0(sample_ids,"-pico")]
  colnames(Plasma_forcor) <- gsub("-pico","",colnames(Plasma_forcor))
  
  forcor1 <- Tissue_T_forcor
  forcor2 <- Plasma_forcor
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,0.001)
  
  forcor1 <- Tissue_N_forcor
  forcor2 <- Plasma_forcor
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/TumorAdjacentNormal_vs_Plasma/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/TumorAdjacentNormal_vs_Plasma/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,0.001)
  
  forcor1 <- (Tissue_T_forcor+1)/(Tissue_N_forcor+1)
  forcor2 <- (Plasma_forcor+1)/(rowSums(Plasma[,-grep("CRC|STAD|NC-PKU-mix.-pico",colnames(Plasma))])+1)
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma_TvsN/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma_TvsN/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,0.001)
}

#Paired correlation between tissue and plasma (pathway level, GSVA)
{
  library(GSVA)
  library(GSA)
  Plasma_genecount <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Plasma_TPM.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  Tissue_genecount <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Tissue_TPM.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  #gset.idx.list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/KEGG_pathway_genes/PATH_ID_NAME_modified.csv")
  gset.idx.list <- GSA::GSA.read.gmt("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/KEGG_pathway_genes/ensembl.gmt")
  names(gset.idx.list$genesets) <- gset.idx.list$geneset.descriptions
  
  Plasma_genecount$Gene <- as.character(lapply(strsplit(rownames(Plasma_genecount),".",fixed = TRUE), function(x) x[1]))
  Plasma_genecount <- aggregate(Plasma_genecount[,-which(colnames(Plasma_genecount)=="Gene")], list(Gene=as.character(Plasma_genecount$Gene)), FUN = sum)
  rownames(Plasma_genecount) <- Plasma_genecount$Gene
  Plasma_genecount <- Plasma_genecount[,-which(colnames(Plasma_genecount)=="Gene")]
  Plasma_genecount <- Plasma_genecount[grep("ENSG",rownames(Plasma_genecount)),]
  
  Tissue_genecount$Gene <- as.character(lapply(strsplit(rownames(Tissue_genecount),".",fixed = TRUE), function(x) x[1]))
  Tissue_genecount <- aggregate(Tissue_genecount[,-which(colnames(Tissue_genecount)=="Gene")], list(Gene=as.character(Tissue_genecount$Gene)), FUN = sum)
  rownames(Tissue_genecount) <- Tissue_genecount$Gene
  Tissue_genecount <- Tissue_genecount[,-which(colnames(Tissue_genecount)=="Gene")]
  Tissue_genecount <- Tissue_genecount[grep("ENSG",rownames(Tissue_genecount)),]
  
  avgexppergene <- rowMeans(Plasma_genecount)
  #plot(density(avgexppergene), xlab="Gene average expression level")
  mask <- avgexppergene > 1
  Plasma_genecount.filt <- Plasma_genecount[mask, ]
  avgexppergene <- rowMeans(Tissue_genecount)
  #plot(density(avgexppergene), xlab="Gene average expression level")
  mask <- avgexppergene > 1
  Tissue_genecount.filt <- Tissue_genecount[mask, ]


  Plasma_pathwayactivity <- gsva(as.matrix(Plasma_genecount.filt), gset.idx.list$genesets, 
       method=c("gsva"),
       kcdf=c("Gaussian"),
       abs.ranking=FALSE,
       min.sz=0,
       max.sz=500,
       parallel.sz=1,
       parallel.type="SOCK",
       mx.diff=TRUE,
       tau=1,
       ssgsea.norm=TRUE,
       verbose=TRUE)
  Tissue_pathwayactivity <- gsva(as.matrix(Tissue_genecount.filt), gset.idx.list$genesets, 
                                 method=c("gsva"),
                                 kcdf=c("Gaussian"),
                                 abs.ranking=FALSE,
                                 min.sz=0,
                                 max.sz=500,
                                 parallel.sz=1,
                                 parallel.type="SOCK",
                                 mx.diff=TRUE,
                                 tau=1,
                                 ssgsea.norm=TRUE,
                                 verbose=TRUE)
  
  Plasma <- Plasma_pathwayactivity 
  Tissue <- Tissue_pathwayactivity
  
  #Plasma <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Plasma_pathwayactivity.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  #Tissue <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Tissue_pathwayactivity.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  
  gene_ids <- rownames(Plasma)
  
  sample_ids <- c("CRC-PKU-27","CRC-PKU-28","CRC-PKU-29","CRC-PKU-30","CRC-PKU-32","CRC-PKU-34",
                  "CRC-PKU-35","CRC-PKU-36","CRC-PKU-37","CRC-PKU-38","CRC-PKU-39","CRC-PKU-40","CRC-PKU-41")
  
  Tissue_T_forcor <- Tissue[gene_ids,paste0(sample_ids,"-T")]
  colnames(Tissue_T_forcor) <- gsub("-T","",colnames(Tissue_T_forcor))
  Tissue_N_forcor <- Tissue[gene_ids,paste0(sample_ids,"-N")]
  colnames(Tissue_N_forcor) <- gsub("-N","",colnames(Tissue_N_forcor))
  Plasma_forcor <- Plasma[gene_ids,paste0(sample_ids,"-pico")]
  colnames(Plasma_forcor) <- gsub("-pico","",colnames(Plasma_forcor))
  
  forcor1 <- Tissue_T_forcor
  forcor2 <- Plasma_forcor
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma_pathway/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma_pathway/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir)
  
  forcor1 <- Tissue_N_forcor
  forcor2 <- Plasma_forcor
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/TumorAdjacentNormal_vs_Plasma_pathway/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/TumorAdjacentNormal_vs_Plasma_pathway/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir)
  
  forcor1 <- (Tissue_T_forcor)-(Tissue_N_forcor)
  forcor2 <- (Plasma_forcor)#-(rowSums(Plasma[,-grep("CRC|STAD|NC-PKU-mix.-pico",colnames(Plasma))]))
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma_TvsN_pathway/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma_TvsN_pathway/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir)
}

#Paired correlation between tissue and plasma (pathway level, mean of pathway genes expression)
{
    #KEGG pathway genes
    PATH_ID_NAME <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/09.Published_data/pathway/PATH_ID_NAME_modified.csv",header = TRUE,row.names = 1)
    #multiomics plasma
    counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Plasma_TPM.txt",sep="\t",header = T, row.names = 1)
    #multiomics tissue
    counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Tissue_TPM.txt",sep="\t",header = T, row.names = 1)
  
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
    write.table(pathway_count,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Tissue_pathway_mean.txt",sep = "\t", quote = FALSE)

    
    Plasma <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Plasma_pathway_mean.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
    Tissue <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Tissue_pathway_mean.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
    
    gene_ids <- rownames(Plasma)
    
    sample_ids <- c("CRC-PKU-27","CRC-PKU-28","CRC-PKU-29","CRC-PKU-30","CRC-PKU-32","CRC-PKU-34",
                    "CRC-PKU-35","CRC-PKU-36","CRC-PKU-37","CRC-PKU-38","CRC-PKU-39","CRC-PKU-40","CRC-PKU-41")
    
    Tissue_T_forcor <- Tissue[gene_ids,paste0(sample_ids,"-T")]
    colnames(Tissue_T_forcor) <- gsub("-T","",colnames(Tissue_T_forcor))
    Tissue_N_forcor <- Tissue[gene_ids,paste0(sample_ids,"-N")]
    colnames(Tissue_N_forcor) <- gsub("-N","",colnames(Tissue_N_forcor))
    Plasma_forcor <- Plasma[gene_ids,paste0(sample_ids,"-pico")]
    colnames(Plasma_forcor) <- gsub("-pico","",colnames(Plasma_forcor))
    
    forcor1 <- Tissue_T_forcor
    forcor2 <- Plasma_forcor
    cor_methods <- "spearman"
    output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma_pathway_mean/"
    dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma_pathway_mean/",recursive = TRUE)
    paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,0.1)
    
    forcor1 <- Tissue_N_forcor
    forcor2 <- Plasma_forcor
    cor_methods <- "spearman"
    output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/TumorAdjacentNormal_vs_Plasma_pathway_mean/"
    dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/TumorAdjacentNormal_vs_Plasma_pathway_mean/",recursive = TRUE)
    paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,0.1)
    
    forcor1 <- (Tissue_T_forcor+1)/(Tissue_N_forcor+1)
    forcor2 <- (Plasma_forcor+1)/(rowSums(Plasma[,-grep("CRC|STAD|NC-PKU-mix.-pico",colnames(Plasma))])+1)
    cor_methods <- "spearman"
    output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma_TvsN_pathway_mean/"
    dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Tumor_vs_Plasma_TvsN_pathway_mean/",recursive = TRUE)
    paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,0.1)
}

#Paired correlation between tissue and plasma (splicing site level)
{
  Plasma <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Splicing/multiomics_paired_CRC_vs_Healthy_Inclevel_matrix.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  Tissue <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Splicing/multiomics_tissue_Tumor_vs_Normal_Inclevel_matrix.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  colnames(Plasma) <- gsub(".","-",fixed = TRUE,colnames(Plasma))
  colnames(Tissue) <- gsub(".","-",fixed = TRUE,colnames(Tissue))
  
  Differential_plasma <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Splicing/multiomics_paired_CRC_vs_Healthy.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  Differential_tissue <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Splicing/multiomics_tissue_Tumor_vs_Normal.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
    
  
  gene_ids <- intersect(rownames(Differential_tissue[Differential_tissue$FDR<0.1,]),rownames(Differential_plasma[Differential_plasma$FDR<0.1,]))
  
  sample_ids <- c("CRC-PKU-27","CRC-PKU-28","CRC-PKU-29","CRC-PKU-30","CRC-PKU-32","CRC-PKU-34",
                  "CRC-PKU-35","CRC-PKU-36","CRC-PKU-37","CRC-PKU-38","CRC-PKU-39","CRC-PKU-40","CRC-PKU-41")
  
  Tissue_T_forcor <- Tissue[gene_ids,paste0(sample_ids,"-T")]
  colnames(Tissue_T_forcor) <- gsub("-T","",colnames(Tissue_T_forcor))
  Tissue_N_forcor <- Tissue[gene_ids,paste0(sample_ids,"-N")]
  colnames(Tissue_N_forcor) <- gsub("-N","",colnames(Tissue_N_forcor))
  Plasma_forcor <- Plasma[gene_ids,paste0(sample_ids,"-pico")]
  colnames(Plasma_forcor) <- gsub("-pico","",colnames(Plasma_forcor))
  
  forcor1 <- Tissue_T_forcor
  forcor2 <- Plasma_forcor
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Splicing_Tumor_vs_Plasma/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Splicing_Tumor_vs_Plasma/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,0.1)
  
  forcor1 <- Tissue_N_forcor
  forcor2 <- Plasma_forcor
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Splicing_TumorAdjacentNormal_vs_Plasma/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Splicing_TumorAdjacentNormal_vs_Plasma/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,0.1)
  
  forcor1 <- (Tissue_T_forcor)-(Tissue_N_forcor)
  forcor2 <- (Plasma_forcor)
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Splicing_Tumor_vs_Plasma_TvsN/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Splicing_Tumor_vs_Plasma_TvsN/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,0.1)
}

#Paired correlation between tissue and plasma (chimric RNA level)
{
  Plasma <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Chimric_RNA/Plasma_JunctionReadCount.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  Tissue <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Chimric_RNA/Tissue_JunctionReadCount.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)

  gene_ids <- intersect(rownames(Plasma),rownames(Tissue))
  
  sample_ids <- c("CRC-PKU-27","CRC-PKU-28","CRC-PKU-29","CRC-PKU-30","CRC-PKU-32","CRC-PKU-34",
                  "CRC-PKU-35","CRC-PKU-36","CRC-PKU-37","CRC-PKU-38","CRC-PKU-39","CRC-PKU-40","CRC-PKU-41")
  
  Tissue_T_forcor <- Tissue[gene_ids,paste0(sample_ids,"-T")]
  colnames(Tissue_T_forcor) <- gsub("-T","",colnames(Tissue_T_forcor))
  Tissue_N_forcor <- Tissue[gene_ids,paste0(sample_ids,"-N")]
  colnames(Tissue_N_forcor) <- gsub("-N","",colnames(Tissue_N_forcor))
  Plasma_forcor <- Plasma[gene_ids,paste0(sample_ids,"-pico")]
  colnames(Plasma_forcor) <- gsub("-pico","",colnames(Plasma_forcor))
  
  forcor1 <- Tissue_T_forcor
  forcor2 <- Plasma_forcor
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Chimeric_Tumor_vs_Plasma/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Chimeric_Tumor_vs_Plasma/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,1)
  
  forcor1 <- Tissue_N_forcor
  forcor2 <- Plasma_forcor
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Chimeric_TumorAdjacentNormal_vs_Plasma/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Chimeric_TumorAdjacentNormal_vs_Plasma/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,1)
  
  forcor1 <- (Tissue_T_forcor)-(Tissue_N_forcor)
  forcor2 <- (Plasma_forcor)
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Chimeric_Tumor_vs_Plasma_TvsN/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Chimeric_Tumor_vs_Plasma_TvsN/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,1)
}

#Paired correlation between tissue and plasma (RNA mutation level)
{
  Plasma <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Mutation/Plasma_Mutation_burden_gene.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  Tissue <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Mutation/Tissue_Mutation_burden_gene.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
  
  gene_ids <- intersect(rownames(Plasma),rownames(Tissue))
  
  sample_ids <- c("CRC-PKU-27","CRC-PKU-28","CRC-PKU-29","CRC-PKU-30","CRC-PKU-32","CRC-PKU-34",
                  "CRC-PKU-35","CRC-PKU-36","CRC-PKU-37","CRC-PKU-38","CRC-PKU-39","CRC-PKU-40","CRC-PKU-41")
  
  Tissue_T_forcor <- Tissue[gene_ids,paste0(sample_ids,"-T")]
  colnames(Tissue_T_forcor) <- gsub("-T","",colnames(Tissue_T_forcor))
  Tissue_N_forcor <- Tissue[gene_ids,paste0(sample_ids,"-N")]
  colnames(Tissue_N_forcor) <- gsub("-N","",colnames(Tissue_N_forcor))
  Plasma_forcor <- Plasma[gene_ids,paste0(sample_ids,"-pico")]
  colnames(Plasma_forcor) <- gsub("-pico","",colnames(Plasma_forcor))
  
  forcor1 <- Tissue_T_forcor
  forcor2 <- Plasma_forcor
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Mutation_Tumor_vs_Plasma/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Mutation_Tumor_vs_Plasma/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,1)
  
  forcor1 <- Tissue_N_forcor
  forcor2 <- Plasma_forcor
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Mutation_TumorAdjacentNormal_vs_Plasma/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Mutation_TumorAdjacentNormal_vs_Plasma/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,1)
  
  forcor1 <- (Tissue_T_forcor)-(Tissue_N_forcor)
  forcor2 <- (Plasma_forcor)
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Mutation_Tumor_vs_Plasma_TvsN/"
  dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Spearman_correlation/Mutation_Tumor_vs_Plasma_TvsN/",recursive = TRUE)
  paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,1)
}

############################################################################
#DNA clinical correlation
{
  #clinical
  {
    clinical <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Sample_clinical_information.csv",header = TRUE,check.names = FALSE,row.names = 1)
    clinical$ID <-rownames(clinical)
  }
  
  #DNA
  {
    mat <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/DNA/CPM-TMM_matrix_gene.correctGC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    mat <- as.data.frame(t(mat))
    mat$ID <- rownames(mat)
    mat$ID <- gsub("-wgs","",mat$ID)
    mat <- mat[-grep("NC",mat$ID),]
  }
  
  forCorrelation <- left_join(clinical,mat,by=c("ID"="ID"))
  rownames(forCorrelation) <- forCorrelation$ID
  forCorrelation <- forCorrelation[,-which(colnames(forCorrelation)=="ID")]
  
  correlation_method <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/CNV_correlation/"
  enmuerate_correlation(forCorrelation,correlation_method,output_dir,ncol(clinical))
}



