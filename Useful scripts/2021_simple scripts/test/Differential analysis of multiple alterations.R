#library preparation
{
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
}

#Function
{
  #function for RNA expression: edgeR exact test
  edgeR_exact_test <- function(mat_raw,des,output_matrix,output_res){
    samples <- des$samples
    #samples <- gsub(".","-",samples,fixed=TRUE)
    #samples <- gsub("X","",samples,fixed=TRUE)
    group <- des$group
    
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
    write.table(mat,output_matrix, sep='\t', quote=FALSE, row.names=TRUE) 
    
    y <- DGEList(counts=mat, samples=samples, group=group)     #规定输入格式，必须包含这几个参数
    keep <- filterByExpr(y,group = group)
    y <- y[keep,]
    y <- calcNormFactors(y, method="TMM")                      #归一化处理 TMM for RNA, CNV and Methylation
    y <- estimateDisp(y)
    test <- exactTest(y, pair = c("negative","positive"), dispersion = "auto")
    #https://support.bioconductor.org/p/64807/
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
    write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE) #输出文件，规定名字
  }
  
  #function for RNA mutation: wilcox rank sum test
  AF_wilcox_test <- function(mat_raw,des,output_matrix,output_res){
    norm_method <- 'NA'
    samples <- des$samples
    #samples <- gsub(".","-",samples,fixed=TRUE)
    #samples <- gsub("X","",samples,fixed=TRUE)
    group <- des$group
    positive <- des[which(des$group=="positive"),]$samples
    negative <- des[which(des$group=="negative"),]$samples
    
    if(norm_method == 'NA' ){
      message('Matrix output without normalization.')
      mat_raw[is.na(mat_raw)] <- 0
      matrix <- mat_raw
    }else{
      message('Matrix output normalized by:',norm_method)
      matrix <- cpm(mat_raw, method=norm_method)
    }
    
    colnames(matrix) <- gsub(".","-",fixed = TRUE,colnames(matrix))
    matrix <- matrix[,as.character(samples)]
    
    test_func <- function(x){
      positive_mean<-mean(x[group=="positive"])
      negative_mean<-mean(x[group=="negative"])
      if (positive_mean == negative_mean) {
        wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='two.sided')$p.value
      } else if (positive_mean > negative_mean) {
        wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='greater')$p.value
      } else if (positive_mean < negative_mean) {
        wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='less')$p.value
      }
    }
    
    #filter Editing events < 20% samples in minial group
    positive_prop <- rowSums(matrix[,positive] > 0)/length(positive)
    negative_prop <- rowSums(matrix[,negative] > 0)/length(negative)
    matrix <- matrix[(negative_prop>=0.2)&(positive_prop>=0.2),]
    write.table(matrix, output_matrix, sep='\t', quote=FALSE, row.names=TRUE)
    
    pvalues <- apply(matrix, 1, test_func)
    matrix_logcpm = log2(matrix + 1)
    logFC <- apply(matrix_logcpm[,positive], 1, mean) -
      apply(matrix_logcpm[,negative], 1, mean)
    deltaAF <- apply(matrix[,positive], 1, mean) -
      apply(matrix[,negative], 1, mean)
    res <- data.frame(log2FoldChange=logFC,
                      deltaAF=deltaAF,
                      pvalue=pvalues, 
                      padj=p.adjust(pvalues, method='BH'),
                      baseMean=apply(matrix, 1, mean))
    write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE)
  }
  
  #function for RNA Editing: wilcox rank sum test
  Editing_ratio_wilcox_test <- function(mat_raw,des,output_matrix,output_res){
    norm_method <- 'NA'
    samples <- des$samples
    #samples <- gsub(".","-",samples,fixed=TRUE)
    #samples <- gsub("X","",samples,fixed=TRUE)
    group <- des$group
    positive <- des[which(des$group=="positive"),]$samples
    negative <- des[which(des$group=="negative"),]$samples
    
    if(norm_method == 'NA' ){
      message('Matrix output without normalization.')
      mat_raw[is.na(mat_raw)] <- 0
      matrix <- mat_raw
    }else{
      message('Matrix output normalized by:',norm_method)
      matrix <- cpm(mat_raw, method=norm_method)
    }
    
    colnames(matrix) <- gsub(".","-",fixed = TRUE,colnames(matrix))
    matrix <- matrix[,as.character(samples)]
    
    
    
    test_func <- function(x){
      positive_mean<-mean(x[group=="positive"])
      negative_mean<-mean(x[group=="negative"])
      if (positive_mean == negative_mean) {
        wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='two.sided')$p.value
      } else if (positive_mean > negative_mean) {
        wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='greater')$p.value
      } else if (positive_mean < negative_mean) {
        wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='less')$p.value
      }
    }
    
    #filter Editing events < 20% samples in minial group
    positive_prop <- rowSums(matrix[,positive] > 0)/length(positive)
    negative_prop <- rowSums(matrix[,negative] > 0)/length(negative)
    matrix <- matrix[(negative_prop>=0.2)&(positive_prop>=0.2),]
    
    write.table(matrix, output_matrix, sep='\t', quote=FALSE, row.names=TRUE)
    
    pvalues <- apply(matrix, 1, test_func)
    matrix_logcpm = log2(matrix + 1)
    logFC <- apply(matrix_logcpm[,positive], 1, mean) -
      apply(matrix_logcpm[,negative], 1, mean)
    deltaAF <- apply(matrix[,positive], 1, mean) -
      apply(matrix[,negative], 1, mean)
    res <- data.frame(log2FoldChange=logFC,
                      deltaAF=deltaAF,
                      pvalue=pvalues, 
                      padj=p.adjust(pvalues, method='BH'),
                      baseMean=apply(matrix, 1, mean))
    write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE)
  }
  
  #function for enrichment analysis for RNA expression: clusterprofiler
  KEGG_GO_Expression <- function(Differential_result,
                                 output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                 pvalue_cutoff,log2Foldchange_cutoff){
    #all_gene <- row.names(Differential_result)
    #strsplit is highly depend on the colnames of alteration matrix
    #all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
    #background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
    #                    filters = "ensembl_gene_id",
    #                    values=all_gene, mart= mart,useCache = FALSE)
    #background <- background[-which(background$entrezgene_id=="")]
    #background <- background$entrezgene_id
    
    filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange < -log2Foldchange_cutoff),])
    filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
    forenrich_down <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                            filters = "ensembl_gene_id",
                            values=filtered_down, mart= mart,useCache = FALSE)
    if(length(which(forenrich_down$entrezgene_id==""))==0){
      forenrich_down <- forenrich_down
    } else {
      forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
    }
    forenrich_down <- forenrich_down$entrezgene_id
    
    filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange > log2Foldchange_cutoff),])
    filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
    forenrich_up <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                          filters = "ensembl_gene_id",
                          values=filtered_up, mart= mart,useCache = FALSE)
    if(length(which(forenrich_up$entrezgene_id==""))==0){
      forenrich_up <- forenrich_up
    } else {
      forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
    }
    forenrich_up <- forenrich_up$entrezgene_id
      
    #KEGG
    {
      KEGG_res_down <- enrichKEGG(
        forenrich_down,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1,
        use_internal_data = FALSE)
      KEGG_output_down <- KEGG_res_down@result
      KEGG_output_down$GeneEnrichedIn <- "Down regulated"
      
      KEGG_res_up <- enrichKEGG(
        forenrich_up,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1,
        use_internal_data = FALSE)
      KEGG_output_up <- KEGG_res_up@result
      KEGG_output_up$GeneEnrichedIn <- "Up regulated"
      
      KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
      write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
    }
    #GO_BP
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
    }
    #GO_CC
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "CC",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "CC",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
    }
    #GO_MF
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "MF",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "MF",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
    }
  }
  
  #function for enrichment analysis for miRNA: clusterprofiler
  KEGG_GO_miRNA <- function(Differential_result,
                            output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                            pvalue_cutoff,log2Foldchange_cutoff){
    #all_tx <- row.names(Differential_result)
    #strsplit is highly depend on the colnames of alteration matrix
    #all_tx <- as.character(lapply(strsplit(all_tx,"\\."),function(x) x[1]))
    #all_mir <- getBM(attributes=c("ensembl_transcript_id", "mirbase_id"),
    #                         filters = "ensembl_transcript_id",
    #                         values=all_tx, mart= mart,useCache = FALSE)
    #all_mir <- all_mir[-which(all_mir$mirbase_id=="")]
    #background <- get_multimir(mirna = all_mir$mirbase_id, summary = TRUE)
    #background <- unique(background@data$target_entrez)
    
    filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange < -log2Foldchange_cutoff),])
    filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
    down_tx <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
    down_mir <- getBM(attributes=c("ensembl_transcript_id", "mirbase_id"),
                     filters = "ensembl_transcript_id",
                     values=down_tx, mart= mart,useCache = FALSE)
    if(length(which(down_mir$entrezgene_id==""))==0){
      down_mir <- down_mir
    } else {
      down_mir <- down_mir[-which(down_mir$mirbase_id==""),]
    }
    forenrich_down <- get_multimir(mirna = down_mir$mirbase_id, summary = TRUE)
    forenrich_down <- unique(forenrich_down@data$target_entrez)
    
    filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange > log2Foldchange_cutoff),])
    filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
    up_tx <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
    up_mir <- getBM(attributes=c("ensembl_transcript_id", "mirbase_id"),
                      filters = "ensembl_transcript_id",
                      values=up_tx, mart= mart,useCache = FALSE)
    if(length(which(up_mir$entrezgene_id==""))==0){
      up_mir <- up_mir
    } else {
      up_mir <- up_mir[-which(up_mir$mirbase_id==""),]
    }
    forenrich_up <- get_multimir(mirna = up_mir$mirbase_id, summary = TRUE)
    forenrich_up <- unique(forenrich_up@data$target_entrez)
    
    #KEGG
    {
      KEGG_res_down <- enrichKEGG(
        forenrich_down,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1,
        use_internal_data = FALSE)
      KEGG_output_down <- KEGG_res_down@result
      KEGG_output_down$GeneEnrichedIn <- "Down regulated"
      
      KEGG_res_up <- enrichKEGG(
        forenrich_up,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1,
        use_internal_data = FALSE)
      KEGG_output_up <- KEGG_res_up@result
      KEGG_output_up$GeneEnrichedIn <- "Up regulated"
      
      KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
      write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
    }
    #GO_BP
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
    }
    #GO_CC
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "CC",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "CC",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
    }
    #GO_MF
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "MF",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "MF",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
    }
  }
  
  #function for enrichment analysis for Splicing: clusterprofiler
  KEGG_GO_Splicing <- function(Differential_result,
                               output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                               PValue_cutoff,IncLevelDifference_cutoff){
    #all_gene <- row.names(Differential_result)
    #strsplit is highly depend on the colnames of alteration matrix
    #all_gene <- as.character(lapply(strsplit(all_gene,"|",fixed = TRUE),function(x) x[2]))
    #all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
    #background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
    #                    filters = "ensembl_gene_id",
    #                    values=all_gene, mart= mart,useCache = FALSE)
    #background <- background[-which(background$entrezgene_id=="")]
    #background <- background$entrezgene_id
    
    filtered_down <- rownames(Differential_result[(Differential_result$PValue<PValue_cutoff)&(Differential_result$IncLevelDifference < -IncLevelDifference_cutoff),])
    filtered_down <- as.character(lapply(strsplit(filtered_down,"|",fixed = TRUE),function(x) x[2]))
    filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
    forenrich_down <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                            filters = "ensembl_gene_id",
                            values=filtered_down, mart= mart,useCache = FALSE)
    if(length(which(forenrich_down$entrezgene_id==""))==0){
      forenrich_down <- forenrich_down
    } else {
      forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
    }
    forenrich_down <- forenrich_down$entrezgene_id
    
    filtered_up <- rownames(Differential_result[(Differential_result$PValue<PValue_cutoff)&(Differential_result$IncLevelDifference > IncLevelDifference_cutoff),])
    filtered_up <- as.character(lapply(strsplit(filtered_up,"|",fixed = TRUE),function(x) x[2]))
    filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
    forenrich_up <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                          filters = "ensembl_gene_id",
                          values=filtered_up, mart= mart,useCache = FALSE)
    if(length(which(forenrich_up$entrezgene_id==""))==0){
      forenrich_up <- forenrich_up
    } else {
      forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
    }
    forenrich_up <- forenrich_up$entrezgene_id
    
    #KEGG
    {
      KEGG_res_down <- enrichKEGG(
        forenrich_down,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1,
        use_internal_data = FALSE)
      KEGG_output_down <- KEGG_res_down@result
      KEGG_output_down$GeneEnrichedIn <- "Down regulated"
      
      KEGG_res_up <- enrichKEGG(
        forenrich_up,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1,
        use_internal_data = FALSE)
      KEGG_output_up <- KEGG_res_up@result
      KEGG_output_up$GeneEnrichedIn <- "Up regulated"
      
      KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
      write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
    }
    #GO_BP
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
    }
    #GO_CC
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "CC",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "CC",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
    }
    #GO_MF
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "MF",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "MF",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
    }
  }
  
  #function for enrichment analysis for Mutation site: clusterprofiler
  KEGG_GO_Mutation <- function(Differential_result,
                               output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                               pvalue_cutoff,deltaAF_cutoff){
    filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF < -deltaAF_cutoff),])
    chr <- unlist(lapply(strsplit(filtered_down,"_"), function(x) x[1]))
    chr <- gsub("chr","",chr)
    start <- as.integer(unlist(lapply(strsplit(filtered_down,"_"), function(x) x[2])))
    if(length(which(is.na(start)))==0){
      chr <- chr
      start <- start
    } else {
      chr <- chr[-which(is.na(start))]
      start <- start[-which(is.na(start))]
    }
    
    i=1
    forenrich_down <- data.frame()
    while(i<=length(chr)){
      message("Processing: ",i,"/",length(chr))
      attributes <- c("ensembl_gene_id","start_position","end_position","entrezgene_id")
      filters <- c("chromosome_name","start","end")
      values <- list(chromosome=chr[i],start=start[i],end=start[i])
      forenrich_down_tmp <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","entrezgene_id"),
                                  filters=c("chromosome_name","start","end"), 
                                  values=values, mart=mart, useCache = FALSE)
      forenrich_down <- rbind(forenrich_down,forenrich_down_tmp)
      i=i+1
    }
    
    forenrich_down <- na.omit(forenrich_down)
    if(length(which(forenrich_down$entrezgene_id==""))==0){
      forenrich_down <- forenrich_down
    } else {
      forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
    }
    forenrich_down <- forenrich_down$entrezgene_id
    forenrich_down <- unique(forenrich_down)
    
    filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF > deltaAF_cutoff),])
    chr <- unlist(lapply(strsplit(filtered_up,"_"), function(x) x[1]))
    chr <- gsub("chr","",chr)
    start <- as.integer(unlist(lapply(strsplit(filtered_up,"_"), function(x) x[2])))
    if(length(which(is.na(start)))==0){
      chr <- chr
      start <- start
    } else {
      chr <- chr[-which(is.na(start))]
      start <- start[-which(is.na(start))]
    }
    
    i=1
    forenrich_up <- data.frame()
    while(i<=length(chr)){
      message("Processing: ",i,"/",length(chr))
      attributes <- c("ensembl_gene_id","start_position","end_position","entrezgene_id")
      filters <- c("chromosome_name","start","end")
      values <- list(chromosome=chr[i],start=start[i],end=start[i])
      forenrich_up_tmp <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","entrezgene_id"),
                                filters=c("chromosome_name","start","end"), 
                                values=values, mart=mart, useCache = FALSE)
      forenrich_up <- rbind(forenrich_up,forenrich_up_tmp)
      i=i+1
    }
    
    forenrich_up <- na.omit(forenrich_up)
    if(length(which(forenrich_up$entrezgene_id==""))==0){
      forenrich_up <- forenrich_up
    } else {
      forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
    }
    forenrich_up <- forenrich_up$entrezgene_id
    forenrich_up <- unique(forenrich_up)
    
    #KEGG
    {
      KEGG_res_down <- enrichKEGG(
        forenrich_down,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1,
        use_internal_data = FALSE)
      KEGG_output_down <- KEGG_res_down@result
      KEGG_output_down$GeneEnrichedIn <- "Down regulated"
      
      KEGG_res_up <- enrichKEGG(
        forenrich_up,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1,
        use_internal_data = FALSE)
      KEGG_output_up <- KEGG_res_up@result
      KEGG_output_up$GeneEnrichedIn <- "Up regulated"
      
      KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
      write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
    }
    #GO_BP
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
    }
    #GO_CC
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "CC",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "CC",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
    }
    #GO_MF
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "MF",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "MF",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
    }
  }
  
  #function for enrichment analysis for Editing site: clusterprofiler
  KEGG_GO_Editing <- function(Differential_result,
                               output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                               pvalue_cutoff,deltaAF_cutoff){
    filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF < -deltaAF_cutoff),])
    chr <- unlist(lapply(strsplit(filtered_down,"_"), function(x) x[1]))
    chr <- gsub("chr","",chr)
    start <- as.integer(unlist(lapply(strsplit(filtered_down,"_"), function(x) x[2])))
    if(length(which(is.na(start)))==0){
      chr <- chr
      start <- start
    } else {
      chr <- chr[-which(is.na(start))]
      start <- start[-which(is.na(start))]
    }
    
    i=1
    forenrich_down <- data.frame()
    while(i<=length(chr)){
      message("Processing: ",i,"/",length(chr))
      attributes <- c("ensembl_gene_id","start_position","end_position","entrezgene_id")
      filters <- c("chromosome_name","start","end")
      values <- list(chromosome=chr[i],start=start[i],end=start[i])
      forenrich_down_tmp <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","entrezgene_id"),
                                  filters=c("chromosome_name","start","end"), 
                                  values=values, mart=mart, useCache = FALSE)
      forenrich_down <- rbind(forenrich_down,forenrich_down_tmp)
      i=i+1
    }
    
    forenrich_down <- na.omit(forenrich_down)
    if(length(which(forenrich_down$entrezgene_id==""))==0){
      forenrich_down <- forenrich_down
    } else {
      forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
    }
    forenrich_down <- forenrich_down$entrezgene_id
    forenrich_down <- unique(forenrich_down)
    
    filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF > deltaAF_cutoff),])
    chr <- unlist(lapply(strsplit(filtered_up,"_"), function(x) x[1]))
    chr <- gsub("chr","",chr)
    start <- as.integer(unlist(lapply(strsplit(filtered_up,"_"), function(x) x[2])))
    if(length(which(is.na(start)))==0){
      chr <- chr
      start <- start
    } else {
      chr <- chr[-which(is.na(start))]
      start <- start[-which(is.na(start))]
    }
    
    i=1
    forenrich_up <- data.frame()
    while(i<=length(chr)){
      message("Processing: ",i,"/",length(chr))
      attributes <- c("ensembl_gene_id","start_position","end_position","entrezgene_id")
      filters <- c("chromosome_name","start","end")
      values <- list(chromosome=chr[i],start=start[i],end=start[i])
      forenrich_up_tmp <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","entrezgene_id"),
                                filters=c("chromosome_name","start","end"), 
                                values=values, mart=mart, useCache = FALSE)
      forenrich_up <- rbind(forenrich_up,forenrich_up_tmp)
      i=i+1
    }
    
    forenrich_up <- na.omit(forenrich_up)
    if(length(which(forenrich_up$entrezgene_id==""))==0){
      forenrich_up <- forenrich_up
    } else {
      forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
    }
    forenrich_up <- forenrich_up$entrezgene_id
    forenrich_up <- unique(forenrich_up)
    
    #KEGG
    {
      KEGG_res_down <- enrichKEGG(
        forenrich_down,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1,
        use_internal_data = FALSE)
      KEGG_output_down <- KEGG_res_down@result
      KEGG_output_down$GeneEnrichedIn <- "Down regulated"
      
      KEGG_res_up <- enrichKEGG(
        forenrich_up,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1,
        use_internal_data = FALSE)
      KEGG_output_up <- KEGG_res_up@result
      KEGG_output_up$GeneEnrichedIn <- "Up regulated"
      
      KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
      write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
    }
    #GO_BP
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
    }
    #GO_CC
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "CC",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "CC",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
    }
    #GO_MF
    {
      library(clusterProfiler)
      GO_res_down <- enrichGO(
        forenrich_down,
        'org.Hs.eg.db',
        ont = "MF",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_down <- GO_res_down@result
      GO_output_down$GeneEnrichedIn <- "Down regulated"
      
      GO_res_up <- enrichGO(
        forenrich_up,
        'org.Hs.eg.db',
        ont = "MF",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        #universe=as.character(background),
        minGSSize = 0,
        maxGSSize = 500,
        qvalueCutoff = 1)
      GO_output_up <- GO_res_up@result
      GO_output_up$GeneEnrichedIn <- "Up regulated"
      
      GO_output <- rbind(GO_output_up,GO_output_down)
      
      write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
    }
  }
}

#Differential analysis and enrichment analysis
{
  ##Expression
  #RNA expression matrix: /data/taoyuhuan/projects/exOmics_RNA/multiomics_paired/output/Intron-spanning/count/count_matrix/featurecount_intron_spanning.txt
  #Differential expression 
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
    #Differential Analysis: gastrointestinal cancer vs healthy donors
    {
    mat_raw <- read.table("./Expression/matrix/featurecount_intron_spanning.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
    mat_raw <- mat_raw[grep("ENSG",rownames(mat_raw)),]
    des <- read.csv("./Expression/group/des_Expression_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
    output_matrix <- "./Expression/matrix/Expression_for_GIvsNC.txt"
    output_res <- "./Expression/output/Expression_GIvsNC_edger_exact.txt"
    edgeR_exact_test(mat_raw,des,output_matrix,output_res)
    }
    #Differential Analysis: STAD vs CRC
    {
      mat_raw <- read.table("./Expression/matrix/featurecount_intron_spanning.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
      mat_raw <- mat_raw[grep("ENSG",rownames(mat_raw)),]
      des <- read.csv("./Expression/group/des_Expression_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./Expression/matrix/Expression_for_STADvsCRC.txt"
      output_res <- "./Expression/output/Expression_STADvsCRC_edger_exact.txt"
      edgeR_exact_test(mat_raw,des,output_matrix,output_res)
    }
    #Differential Analysis: STAD vs NC
    {
      mat_raw <- read.table("./Expression/matrix/featurecount_intron_spanning.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
      mat_raw <- mat_raw[grep("ENSG",rownames(mat_raw)),]
      des <- read.csv("./Expression/group/des_Expression_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./Expression/matrix/Expression_for_STADvsNC.txt"
      output_res <- "./Expression/output/Expression_STADvsNC_edger_exact.txt"
      edgeR_exact_test(mat_raw,des,output_matrix,output_res)
    }
    #Differential Analysis: CRC vs NC
    {
      mat_raw <- read.table("./Expression/matrix/featurecount_intron_spanning.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
      mat_raw <- mat_raw[grep("ENSG",rownames(mat_raw)),]
      des <- read.csv("./Expression/group/des_Expression_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./Expression/matrix/Expression_for_CRCvsNC.txt"
      output_res <- "./Expression/output/Expression_CRCvsNC_edger_exact.txt"
      edgeR_exact_test(mat_raw,des,output_matrix,output_res)
    }
  }
  #Enrichment analysis
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
    #library biomart
    {
     library(biomaRt)
     mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    }
    #Enrichment analysis: gastrointestinal cancer vs healthy donors
    {
    Differential_result <- read.csv("./Expression/output/Expression_GIvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
    output_KEGG <- "./Expression/output/Expression_GIvsNC_KEGG.txt"
    output_GO_BP <- "./Expression/output/Expression_GIvsNC_GO_BP.txt"
    output_GO_CC <- "./Expression/output/Expression_GIvsNC_GO_CC.txt"
    output_GO_MF <- "./Expression/output/Expression_GIvsNC_GO_MF.txt"
    pvalue_cutoff <- 0.05
    log2Foldchange_cutoff <- 0.59
  
    KEGG_GO_Expression(Differential_result,
                      output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                      pvalue_cutoff,log2Foldchange_cutoff)
    }
    #Enrichment analysis: STAD vs CRC
    {
      Differential_result <- read.csv("./Expression/output/Expression_STADvsCRC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Expression/output/Expression_STADvsCRC_KEGG.txt"
      output_GO_BP <- "./Expression/output/Expression_STADvsCRC_GO_BP.txt"
      output_GO_CC <- "./Expression/output/Expression_STADvsCRC_GO_CC.txt"
      output_GO_MF <- "./Expression/output/Expression_STADvsCRC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      log2Foldchange_cutoff <- 0.59
      
      KEGG_GO_Expression(Differential_result,
                         output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                         pvalue_cutoff,log2Foldchange_cutoff)
    }
    #Enrichment analysis: STAD vs NC
    {
      Differential_result <- read.csv("./Expression/output/Expression_STADvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Expression/output/Expression_STADvsNC_KEGG.txt"
      output_GO_BP <- "./Expression/output/Expression_STADvsNC_GO_BP.txt"
      output_GO_CC <- "./Expression/output/Expression_STADvsNC_GO_CC.txt"
      output_GO_MF <- "./Expression/output/Expression_STADvsNC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      log2Foldchange_cutoff <- 0.59
      
      KEGG_GO_Expression(Differential_result,
                         output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                         pvalue_cutoff,log2Foldchange_cutoff)
    }
    #Enrichment analysis: CRC vs NC
    {
      Differential_result <- read.csv("./Expression/output/Expression_CRCvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Expression/output/Expression_CRCvsNC_KEGG.txt"
      output_GO_BP <- "./Expression/output/Expression_CRCvsNC_GO_BP.txt"
      output_GO_CC <- "./Expression/output/Expression_CRCvsNC_GO_CC.txt"
      output_GO_MF <- "./Expression/output/Expression_CRCvsNC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      log2Foldchange_cutoff <- 0.59
      
      KEGG_GO_Expression(Differential_result,
                         output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                         pvalue_cutoff,log2Foldchange_cutoff)
    }
  }
  
  ##miRNA
  #miRNA expression matrix: /data/taoyuhuan/projects/small_RNAseq_pipeline/output/multiomics_qia/count_matrix
  #Differential expression
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
    #Differential Analysis: gastrointestinal cancer vs healthy donors
    {
      mat_raw <- read.table("./miRNA/matrix/miRNA.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
      des <- read.csv("./miRNA/group/des_miRNA_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./miRNA/matrix/miRNA_for_GIvsNC.txt"
      output_res <- "./miRNA/output/miRNA_GIvsNC_edger_exact.txt"
      edgeR_exact_test(mat_raw,des,output_matrix,output_res)
    }
    #Differential Analysis: STAD vs CRC
    {
      mat_raw <- read.table("./miRNA/matrix/miRNA.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
      des <- read.csv("./miRNA/group/des_miRNA_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./miRNA/matrix/miRNA_for_STADvsCRC.txt"
      output_res <- "./miRNA/output/miRNA_STADvsCRC_edger_exact.txt"
      edgeR_exact_test(mat_raw,des,output_matrix,output_res)
    }
    #Differential Analysis: STAD vs NC
    {
      mat_raw <- read.table("./miRNA/matrix/miRNA.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
      des <- read.csv("./miRNA/group/des_miRNA_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./miRNA/matrix/miRNA_for_STADvsNC.txt"
      output_res <- "./miRNA/output/miRNA_STADvsNC_edger_exact.txt"
      edgeR_exact_test(mat_raw,des,output_matrix,output_res)
    }
    #Differential Analysis: CRC vs NC
    {
      mat_raw <- read.table("./miRNA/matrix/miRNA.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
      des <- read.csv("./miRNA/group/des_miRNA_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./miRNA/matrix/miRNA_for_CRCvsNC.txt"
      output_res <- "./miRNA/output/miRNA_CRCvsNC_edger_exact.txt"
      edgeR_exact_test(mat_raw,des,output_matrix,output_res)
    }
  }
  
  #Enrichment analysis
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
    #library biomart
    {
      library(biomaRt)
      mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    }
    #Enrichment analysis: gastrointestinal cancer vs healthy donors
    {
      Differential_result <- read.csv("./miRNA/output/miRNA_GIvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./miRNA/output/miRNA_GIvsNC_KEGG.txt"
      output_GO_BP <- "./miRNA/output/miRNA_GIvsNC_GO_BP.txt"
      output_GO_CC <- "./miRNA/output/miRNA_GIvsNC_GO_CC.txt"
      output_GO_MF <- "./miRNA/output/miRNA_GIvsNC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      log2Foldchange_cutoff <- 0.59
      
      KEGG_GO_miRNA(Differential_result,
                         output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                         pvalue_cutoff,log2Foldchange_cutoff)
    }
    #Enrichment analysis: STAD vs CRC
    {
      Differential_result <- read.csv("./miRNA/output/miRNA_STADvsCRC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./miRNA/output/miRNA_STADvsCRC_KEGG.txt"
      output_GO_BP <- "./miRNA/output/miRNA_STADvsCRC_GO_BP.txt"
      output_GO_CC <- "./miRNA/output/miRNA_STADvsCRC_GO_CC.txt"
      output_GO_MF <- "./miRNA/output/miRNA_STADvsCRC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      log2Foldchange_cutoff <- 0.59
      
      KEGG_GO_miRNA(Differential_result,
                         output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                         pvalue_cutoff,log2Foldchange_cutoff)
    }
    #Enrichment analysis: STAD vs NC
    {
      Differential_result <- read.csv("./miRNA/output/miRNA_STADvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./miRNA/output/miRNA_STADvsNC_KEGG.txt"
      output_GO_BP <- "./miRNA/output/miRNA_STADvsNC_GO_BP.txt"
      output_GO_CC <- "./miRNA/output/miRNA_STADvsNC_GO_CC.txt"
      output_GO_MF <- "./miRNA/output/miRNA_STADvsNC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      log2Foldchange_cutoff <- 0.59
      
      KEGG_GO_miRNA(Differential_result,
                    output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                    pvalue_cutoff,log2Foldchange_cutoff)
    }
    #Enrichment analysis: CRC vs NC
    {
      Differential_result <- read.csv("./miRNA/output/miRNA_CRCvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./miRNA/output/miRNA_CRCvsNC_KEGG.txt"
      output_GO_BP <- "./miRNA/output/miRNA_CRCvsNC_GO_BP.txt"
      output_GO_CC <- "./miRNA/output/miRNA_CRCvsNC_GO_CC.txt"
      output_GO_MF <- "./miRNA/output/miRNA_CRCvsNC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      log2Foldchange_cutoff <- 0.59
      
      KEGG_GO_miRNA(Differential_result,
                    output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                    pvalue_cutoff,log2Foldchange_cutoff)
    }
  }
  
  ##Splicing
  #splicing output: /data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/multiomics_paired/allpassed-20211020
  #example inclevel matrix: /data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/multiomics_paired/allpassed-20211020/CRCvsNC/summary/Splicing_CRCvsNC_allpassed.txt
  #Differential splicing done by rMATs
  #example differential splicing: /data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/multiomics_paired/allpassed-20211020/CRCvsNC/summary/Splicing_differential_CRCvsNC_allpassed.txt
  #readin stats and inclevel(for example, Rscript is located in /data/taoyuhuan/projects/exOmics_RNA/level_3_splicing)
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
    file_tosum <- "./Splicing/output"
    files <- dir(file_tosum)
    stats <- grep("stats",files,value = TRUE)
    inclevels <- grep("inc_level",files,value = TRUE)
    output_stat <- "./Splicing/output/Splicing_GIvsNC_rMATs.txt"
    output_inclevel <- "./Splicing/matrix/Inclevel_for_GIvsNC.txt"
    
    stat_all <- data.frame()
    for(i in stats){
      stat <- read.csv(paste0(file_tosum,"/",i),sep = "\t", header = TRUE, row.names = 1)
      stat_all <- rbind(stat_all,stat)
    }
    write.table(stat_all,output_stat,quote = FALSE, sep = "\t")
    
    inclevel_matrix <- data.frame()
    for(i in inclevels){
      inclevel <- read.csv(paste0(file_tosum,"/",i),sep = "\t", header = TRUE, row.names = 1)
      inclevel_matrix <- rbind(inclevel_matrix,inclevel)
    }
    write.table(inclevel_matrix,output_inclevel,quote = FALSE, sep = "\t")
    
  }
  
  #Enrichment analysis
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
    #library biomart
    {
      library(biomaRt)
      mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    }
    #Enrichment analysis: gastrointestinal cancer vs healthy donors
    {
      Differential_result <- read.csv("./Splicing/output/Splicing_differential_GIvsNC_allpassed.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Splicing/output/Splicing_GIvsNC_KEGG.txt"
      output_GO_BP <- "./Splicing/output/Splicing_GIvsNC_GO_BP.txt"
      output_GO_CC <- "./Splicing/output/Splicing_GIvsNC_GO_CC.txt"
      output_GO_MF <- "./Splicing/output/Splicing_GIvsNC_GO_MF.txt"
      PValue_cutoff <- 0.05
      IncLevelDifference_cutoff <- 0.05
      
      KEGG_GO_Splicing(Differential_result,
                    output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                    PValue_cutoff,IncLevelDifference_cutoff)
    }
    #Enrichment analysis: CRC vs NC
    {
      Differential_result <- read.csv("./Splicing/output/Splicing_differential_CRCvsNC_allpassed.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Splicing/output/Splicing_CRCvsNC_KEGG.txt"
      output_GO_BP <- "./Splicing/output/Splicing_CRCvsNC_GO_BP.txt"
      output_GO_CC <- "./Splicing/output/Splicing_CRCvsNC_GO_CC.txt"
      output_GO_MF <- "./Splicing/output/Splicing_CRCvsNC_GO_MF.txt"
      PValue_cutoff <- 0.05
      IncLevelDifference_cutoff <- 0.05
      
      KEGG_GO_Splicing(Differential_result,
                       output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                       PValue_cutoff,IncLevelDifference_cutoff)
    }
    #Enrichment analysis: STAD vs NC
    {
      Differential_result <- read.csv("./Splicing/output/Splicing_differential_STADvsNC_allpassed.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Splicing/output/Splicing_STADvsNC_KEGG.txt"
      output_GO_BP <- "./Splicing/output/Splicing_STADvsNC_GO_BP.txt"
      output_GO_CC <- "./Splicing/output/Splicing_STADvsNC_GO_CC.txt"
      output_GO_MF <- "./Splicing/output/Splicing_STADvsNC_GO_MF.txt"
      PValue_cutoff <- 0.05
      IncLevelDifference_cutoff <- 0.05
      
      KEGG_GO_Splicing(Differential_result,
                       output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                       PValue_cutoff,IncLevelDifference_cutoff)
    }
    #Enrichment analysis: STAD vs CRC
    {
      Differential_result <- read.csv("./Splicing/output/Splicing_differential_STADvsCRC_allpassed.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Splicing/output/Splicing_STADvsCRC_KEGG.txt"
      output_GO_BP <- "./Splicing/output/Splicing_STADvsCRC_GO_BP.txt"
      output_GO_CC <- "./Splicing/output/Splicing_STADvsCRC_GO_CC.txt"
      output_GO_MF <- "./Splicing/output/Splicing_STADvsCRC_GO_MF.txt"
      PValue_cutoff <- 0.05
      IncLevelDifference_cutoff <- 0.05
      
      KEGG_GO_Splicing(Differential_result,
                       output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                       PValue_cutoff,IncLevelDifference_cutoff)
    }
  }
  
  ##APA(too less differential APA)
  #APA output: /data/taoyuhuan/projects/exOmics_RNA/level_3_APA/multiomics_paired/matrix
  #Enrichment analysis
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
  }
  
  ##mutation(gene based uEMD method to discover differential mutated genes does not enrich cancer genes)
  #mutation output (mutation burden per gene): /data/taoyuhuan/projects/Machine_learning/20210722_multiomics/SNP_20211011/Matrix/Mutation_burden_gene.txt
  #Differential mutation analysis by uEMDscore
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
    mutation_matrix <- read.csv("./Mutation/Mutation_burden_gene-allpassed.txt",sep = "\t",header = TRUE, row.names = 1)
    
    #rank normalize mutaton or variation counts
    fastRank = function(x) { 
      x[x!=0] = rank(x[x!=0], ties.method="min")+length(which(x==0)) 
      x/length(x) 
    }
    
    #bin counts to build histogram
    bins = function(v, p=100) {
      l = length(v)
      nBins = rep(0,p)
      for(val in v){
        nBins[ceiling(val*(p-1))+1] = nBins[ceiling(val*(p-1))+1]+1
      }
      nBins = nBins/l
      nBins
    }
    
    #compute unidirectional EMD between mutationas and variations
    uEMD = function(tBins, nBins, p=100) {
      sum = 0
      move = 0
      for(i in p:2){
        move = move+tBins[i]-nBins[i]
        sum = sum+max(move,0)
      }
      sum
    }
    
    #compute unidirectional EMD between mutationas and variations
    uEMD = function(tBins, nBins, p=100) {
      sum = 0
      move = 0
      for(i in p:2){
        move = move+tBins[i]-nBins[i]
        sum = sum+max(move,0)
      }
      sum
    }
    
    #generate random uEMDs to compute FDRs
    generateRandomEMDs = function(tRank, nRankBinned) {
      permRank = t(apply(tRank, 1, sample))
      permRankBinned = apply(permRank, 2, bins)
      randEMDs = sapply(1:dim(nRankBinned)[2], function(x) uEMD(permRankBinned[,x], nRankBinned[,x]))
      randEMDs
    }
    
    #compute FDRs based on random uEMDs
    computeFDR = function(uEMDscores, randEMDs) {
      FDRs = sapply(uEMDscores, function(x) length(which(randEMDs>=x))/(length(which(uEMDscores>=x))))
      FDRs
    }
    
    #compute q-values from FDRs
    computeQ = function(FDRs, uEMDscores) {
      qVal = sapply(1:length(FDRs), function(x) min(FDRs[uEMDscores<=uEMDscores[x]]))
      qVal
    }
    
    
    rank_matrix <- as.data.frame(apply(mutation_matrix, MARGIN = 2, FUN = fastRank))
    
    cancer_rank <- as.data.frame(t(rank_matrix[,grep("CRC|STAD",colnames(rank_matrix))]))
    healthy_rank <- as.data.frame(t(rank_matrix[,grep("NC",colnames(rank_matrix))]))
    
    cancer_bins <- lapply(cancer_rank, FUN = bins)
    healthy_bins <- lapply(healthy_rank, FUN = bins)
    
    uEMDscore = sapply(1:length(cancer_bins), 
                       function(x) uEMD(as.numeric(unlist(cancer_bins[x])), as.numeric(unlist(healthy_bins[x]))))
    p=5
    #compute q-values, p determines number of times random uEMDs are generated
    FDRs = rowSums(sapply(1:p, function(x) 
      computeFDR(uEMDscore, generateRandomEMDs(as.matrix(cancer_rank), as.data.frame(healthy_bins)))))/p
    qVal = computeQ(FDRs, uEMDscore)
    
    
    diff_mutation <- data.frame(names(cancer_bins),uEMDscore,FDRs,qVal)
    
    violin_cancer <- data.frame(distribution = cancer_bins$`ENSG00000089009|RPL6`, group = "cancer")
    violin_healthy <- data.frame(distribution = healthy_bins$`ENSG00000089009|RPL6`, group = "healthy")
    
    violin <- rbind(violin_cancer,violin_healthy)
    
    ggplot(violin,aes(x=group,y=distribution)) + geom_violin()
  }
  
  #mutation output (variant allele fraction per site):/data/taoyuhuan/projects/Machine_learning/20210722_multiomics/SNP_20211011/Matrix/final_SNP_AF.txt
  #Differential mutation analysis by allele fraction per site by wilcox rank sum test
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
    #Differential Analysis: gastrointestinal cancer vs healthy donors
    {
      mat_raw <- read.csv("./Mutation/matrix/final_SNP_AF.txt",sep = "\t",header = TRUE, row.names = 1)
      des <- read.csv("./Mutation/group/des_Mutation_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./Mutation/matrix/Mutation_site_for_GIvsNC.txt"
      output_res <- "./Mutation/output/Mutation_site_GIvsNC_wilcox.txt"
      AF_wilcox_test(mat_raw,des,output_matrix,output_res)
    }
    #Differential Analysis: CRC vs NC
    {
      mat_raw <- read.csv("./Mutation/matrix/final_SNP_AF.txt",sep = "\t",header = TRUE, row.names = 1)
      des <- read.csv("./Mutation/group/des_Mutation_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./Mutation/matrix/Mutation_site_for_CRCvsNC.txt"
      output_res <- "./Mutation/output/Mutation_site_CRCvsNC_wilcox.txt"
      AF_wilcox_test(mat_raw,des,output_matrix,output_res)
    }
    #Differential Analysis: STAD vs NC
    {
      mat_raw <- read.csv("./Mutation/matrix/final_SNP_AF.txt",sep = "\t",header = TRUE, row.names = 1)
      des <- read.csv("./Mutation/group/des_Mutation_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./Mutation/matrix/Mutation_site_for_STADvsNC.txt"
      output_res <- "./Mutation/output/Mutation_site_STADvsNC_wilcox.txt"
      AF_wilcox_test(mat_raw,des,output_matrix,output_res)
    }
    #Differential Analysis: STAD vs CRC
    {
      mat_raw <- read.csv("./Mutation/matrix/final_SNP_AF.txt",sep = "\t",header = TRUE, row.names = 1)
      des <- read.csv("./Mutation/group/des_Mutation_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./Mutation/matrix/Mutation_site_for_STADvsCRC.txt"
      output_res <- "./Mutation/output/Mutation_site_STADvsCRC_wilcox.txt"
      AF_wilcox_test(mat_raw,des,output_matrix,output_res)
    }
  }
  
  #Enrichment analysis
  {  
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
    #library biomart
    {
      library(biomaRt)
      mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    }
    
    #Enrichment analysis: gastrointestinal cancer vs healthy donors
    {
      Differential_result <- read.csv("./Mutation/output/Mutation_site_GIvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Mutation/output/Mutation_GIvsNC_KEGG.txt"
      output_GO_BP <- "./Mutation/output/Mutation_GIvsNC_GO_BP.txt"
      output_GO_CC <- "./Mutation/output/Mutation_GIvsNC_GO_CC.txt"
      output_GO_MF <- "./Mutation/output/Mutation_GIvsNC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      deltaAF_cutoff <- 0.2
      
      KEGG_GO_Mutation(Differential_result,
                       output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                       pvalue_cutoff,deltaAF_cutoff)
    }
    #Enrichment analysis: CRC vs NC
    {
      Differential_result <- read.csv("./Mutation/output/Mutation_site_CRCvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Mutation/output/Mutation_CRCvsNC_KEGG.txt"
      output_GO_BP <- "./Mutation/output/Mutation_CRCvsNC_GO_BP.txt"
      output_GO_CC <- "./Mutation/output/Mutation_CRCvsNC_GO_CC.txt"
      output_GO_MF <- "./Mutation/output/Mutation_CRCvsNC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      deltaAF_cutoff <- 0.2
      
      KEGG_GO_Mutation(Differential_result,
                       output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                       pvalue_cutoff,deltaAF_cutoff)
    }
    #Enrichment analysis: STAD vs NC
    {
      Differential_result <- read.csv("./Mutation/output/Mutation_site_STADvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Mutation/output/Mutation_STADvsNC_KEGG.txt"
      output_GO_BP <- "./Mutation/output/Mutation_STADvsNC_GO_BP.txt"
      output_GO_CC <- "./Mutation/output/Mutation_STADvsNC_GO_CC.txt"
      output_GO_MF <- "./Mutation/output/Mutation_STADvsNC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      deltaAF_cutoff <- 0.2
      
      KEGG_GO_Mutation(Differential_result,
                       output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                       pvalue_cutoff,deltaAF_cutoff)
    }
    #Enrichment analysis: STAD vs CRC
    {
      Differential_result <- read.csv("./Mutation/output/Mutation_site_STADvsCRC_wilcox.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Mutation/output/Mutation_STADvsCRC_KEGG.txt"
      output_GO_BP <- "./Mutation/output/Mutation_STADvsCRC_GO_BP.txt"
      output_GO_CC <- "./Mutation/output/Mutation_STADvsCRC_GO_CC.txt"
      output_GO_MF <- "./Mutation/output/Mutation_STADvsCRC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      deltaAF_cutoff <- 0.2
      
      KEGG_GO_Mutation(Differential_result,
                       output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                       pvalue_cutoff,deltaAF_cutoff)
    }
  }  
    
  ##Editing
  #Editing output: /data/taoyuhuan/projects/Machine_learning/20210722_multiomics/EDIT/editing_all_noNA.txt (this contains more than multi omics samples)
  #Differential editing analysis by editing ratio per site by wilcox rank sum test
  {
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
    #Differential Analysis: gastrointestinal cancer vs healthy donors
    {
      mat_raw <- read.csv("./Editing/matrix/editing_all_noNA.txt",sep = "\t",header = TRUE, row.names = 1)
      des <- read.csv("./Editing/group/des_Editing_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./Editing/matrix/Editing_site_for_GIvsNC.txt"
      output_res <- "./Editing/output/Editing_site_GIvsNC_wilcox.txt"
      Editing_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
    }
    
    #Differential Analysis: CRC vs NC
    {
      mat_raw <- read.csv("./Editing/matrix/editing_all_noNA.txt",sep = "\t",header = TRUE, row.names = 1)
      des <- read.csv("./Editing/group/des_Editing_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./Editing/matrix/Editing_site_for_CRCvsNC.txt"
      output_res <- "./Editing/output/Editing_site_CRCvsNC_wilcox.txt"
      Editing_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
    }
    
    #Differential Analysis: STAD vs NC
    {
      mat_raw <- read.csv("./Editing/matrix/editing_all_noNA.txt",sep = "\t",header = TRUE, row.names = 1)
      des <- read.csv("./Editing/group/des_Editing_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./Editing/matrix/Editing_site_for_STADvsNC.txt"
      output_res <- "./Editing/output/Editing_site_STADvsNC_wilcox.txt"
      Editing_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
    }
    
    #Differential Analysis: STAD vs CRC
    {
      mat_raw <- read.csv("./Editing/matrix/editing_all_noNA.txt",sep = "\t",header = TRUE, row.names = 1)
      des <- read.csv("./Editing/group/des_Editing_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
      output_matrix <- "./Editing/matrix/Editing_site_for_STADvsCRC.txt"
      output_res <- "./Editing/output/Editing_site_STADvsCRC_wilcox.txt"
      Editing_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
    }
  }
  #Enrichment analysis
  {  
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
    #library biomart
    {
      library(biomaRt)
      mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    }
    
    #Enrichment analysis: gastrointestinal cancer vs healthy donors
    {
      Differential_result <- read.csv("./Editing/output/Editing_site_GIvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Editing/output/Editing_GIvsNC_KEGG.txt"
      output_GO_BP <- "./Editing/output/Editing_GIvsNC_GO_BP.txt"
      output_GO_CC <- "./Editing/output/Editing_GIvsNC_GO_CC.txt"
      output_GO_MF <- "./Editing/output/Editing_GIvsNC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      deltaAF_cutoff <- 0
      
      KEGG_GO_Editing(Differential_result,
                       output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                       pvalue_cutoff,deltaAF_cutoff)
    }
    #Enrichment analysis: CRC vs NC
    {
      Differential_result <- read.csv("./Editing/output/Editing_site_CRCvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Editing/output/Editing_CRCvsNC_KEGG.txt"
      output_GO_BP <- "./Editing/output/Editing_CRCvsNC_GO_BP.txt"
      output_GO_CC <- "./Editing/output/Editing_CRCvsNC_GO_CC.txt"
      output_GO_MF <- "./Editing/output/Editing_CRCvsNC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      deltaAF_cutoff <- 0
      
      KEGG_GO_Mutation(Differential_result,
                       output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                       pvalue_cutoff,deltaAF_cutoff)
    }
    #Enrichment analysis: STAD vs NC
    {
      Differential_result <- read.csv("./Editing/output/Editing_site_STADvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Editing/output/Editing_STADvsNC_KEGG.txt"
      output_GO_BP <- "./Editing/output/Editing_STADvsNC_GO_BP.txt"
      output_GO_CC <- "./Editing/output/Editing_STADvsNC_GO_CC.txt"
      output_GO_MF <- "./Editing/output/Editing_STADvsNC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      deltaAF_cutoff <- 0
      
      KEGG_GO_Mutation(Differential_result,
                       output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                       pvalue_cutoff,deltaAF_cutoff)
    }
    #Enrichment analysis: STAD vs CRC
    {
      Differential_result <- read.csv("./Editing/output/Editing_site_STADvsCRC_wilcox.txt",header = T, row.names = 1,sep = "\t")
      output_KEGG <- "./Editing/output/Editing_STADvsCRC_KEGG.txt"
      output_GO_BP <- "./Editing/output/Editing_STADvsCRC_GO_BP.txt"
      output_GO_CC <- "./Editing/output/Editing_STADvsCRC_GO_CC.txt"
      output_GO_MF <- "./Editing/output/Editing_STADvsCRC_GO_MF.txt"
      pvalue_cutoff <- 0.05
      deltaAF_cutoff <- 0
      
      KEGG_GO_Mutation(Differential_result,
                       output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                       pvalue_cutoff,deltaAF_cutoff)
    }
  }  
  
  ##ASE
  # ASEP is designed for paired RNA-seq data, not suitable.Paired RNA-seq data is used for two conditions analysis detecting differential ASE, where paired means the same individual is sequenced under both conditions (i.e., pre- vs post-treatment).
  
}
