#enmuerate correlation calculation
{
  ##preparation
  library(pspearman)
  spearman_CI <- function(x, y, alpha = 0.05){
    rs <- cor(x, y, method = "spearman", use = "complete.obs")
    n <- sum(complete.cases(x, y))
    sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
  }
  
  ##initiation
  #output
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation/Plasma and PBMC immune fraction_spearman/"
  #method
  cor_methods <- "spearman"
  #read in data matrix, row: sample, column: numeric score. 
  forcor <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation/Correlation between plasma expression and PBMC immune fraction.csv",header = TRUE, row.names = 1)
  
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = ncol(forcor), clear = FALSE, width= 60)
  
  result_final <- as.data.frame(matrix(numeric(0),ncol=7))
  colnames(result_final) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
  
  result <- as.data.frame(matrix(numeric(0),ncol=7))
  colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
  i=1
  while(i<=(ncol(forcor)-1)){
    #print(i)
    j=1
    while(j<=(ncol(forcor)-i)){
      
      if(cor_methods=="pearson"){
        r_twosided <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods)
        if (is.na(r_twosided$estimate)) {
          r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
        } else if (r_twosided$estimate>0){
          r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "greater")
        } else if (r_twosided$estimate<0) {
          r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "less")
        } else {
          r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
        }
        result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,r$conf.int[1],r$conf.int[2],r_twosided$p.value)
      } else if(cor_methods=="spearman") {
        r_twosided <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
        if (is.na(r_twosided$estimate)) {
          r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
        } else if(r_twosided$estimate>0){
          r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "greater", approximation = "exact")
        } else if (r_twosided$estimate<0) {
          r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "less", approximation = "exact")
        } else {
          r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
        }
        result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,spearman_CI(forcor[,i],forcor[,i+j])[1],spearman_CI(forcor[,i],forcor[,i+j])[2],r_twosided$p.value)
      }
      
      colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval")
      result <- rbind(result,result_tmp)
      j=j+1
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
}