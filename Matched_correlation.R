#matched correlation
{
  ###onfig
  A_forcor <- "matrix A"
  B_forcor <- "matrix B" #matrix A and B should have same colnames and rownames
  output_dir <- "/your/out/put/directory"
  
  ###main function
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = ncol(A_forcor), clear = FALSE, width= 60)
  
  result_final <- as.data.frame(matrix(numeric(0),ncol=7))
  colnames(result_final) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
  
  result <- as.data.frame(matrix(numeric(0),ncol=7))
  colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
  
  i=1
  while(i<=(ncol(A_forcor))){
    r_twosided <- spearman.test(A_forcor[,i],B_forcor[,i],alternative = "two.sided", approximation = "exact")
    if (is.na(r_twosided$estimate)) {
      r <- spearman.test(A_forcor[,i],B_forcor[,i],alternative = "two.sided", approximation = "exact")
    } else if(r_twosided$estimate>0){
      r <- spearman.test(A_forcor[,i],B_forcor[,i],alternative = "greater", approximation = "exact")
    } else if (r_twosided$estimate<0) {
      r <- spearman.test(A_forcor[,i],B_forcor[,i],alternative = "less", approximation = "exact")
    } else {
      r <- spearman.test(A_forcor[,i],B_forcor[,i],alternative = "two.sided", approximation = "exact")
    }
    result_tmp <- data.frame(paste0(colnames(A_forcor)[i],"-plasma vs ",colnames(B_forcor)[i],"-normal"), ######need to manully change suffix, such as "-plasma","-normal"
                             r$estimate,r$p.value,r$method,
                             spearman_CI(A_forcor[,i],B_forcor[,i])[1],spearman_CI(A_forcor[,i],B_forcor[,i])[2],r_twosided$p.value)
    
    colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    result <- rbind(result,result_tmp)
    
    if(i%%5000==0){
      result_final <- rbind(result_final,result)
      write.csv(result,paste0(output_dir,i/5000,"_Result5000.csv"))
      result <- as.data.frame(matrix(numeric(0),ncol=7))
      colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
      i=i+1
    } else if(i==(ncol(A_forcor))){
      result_final <- rbind(result_final,result)
      write.csv(result,paste0(output_dir,i%/%5000+1,"_Result",i%%5000,".csv"))
      i=i+1
    } else {
      i=i+1
    }
    pb$tick()
    Sys.sleep(1 / 100)
  }
  
  result_final$FDR <- p.adjust(result_final$Pvalue,"BH",n=length(result_final$Pvalue))
  
  write.csv(result_final,paste0(output_dir,"Result_final_plasma_vs_Normal.csv"))
  
  result_final <- result_final[grep("ENSG",result_final$DataNames),]
  
  result_final$FDR <- p.adjust(result_final$Pvalue,"BH",n=length(result_final$Pvalue))
  
  write.csv(result_final,paste0(output_dir,"Result_final_plasma_vs_Normal_ensembl.csv"))
  
}