#### environment preparation
{
  library(progress)
  library(dplyr)
  library(pspearman)
}

#functions
{
  ##preparation
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
      
      if(i%%1==0){
        result_final <- rbind(result_final,result)
        write.csv(result,paste0(output_dir,i/1,"_Result20.csv"))
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
    
    #result_final2 <- as.data.frame(matrix(numeric(0),ncol = ncol(forcor),nrow = ))
  }
}

#check correlation
{
  RP_mRNA <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/RP_mRNA_matrix/GSE133684_count_matrix/GSE133684_count_matrix_RP_mRNA.txt",sep = "\t",header = TRUE,row.names = 1)
  TPM <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/GSE133684_count_matrix_TPM.txt",sep = "\t",header = TRUE,row.names = 1)
  
  forcor <- cbind(RP_mRNA,t(TPM))
  cor_methods <- "spearman"
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2022/02.RP_mRNA/correlation/GSE133684"
  dir.create(output_dir,recursive = TRUE)
  enmuerate_correlation(forcor,cor_methods,output_dir,1)
}