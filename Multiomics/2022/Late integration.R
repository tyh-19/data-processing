library(glmnet)
#Late integration of different modals/omics
workdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220629_xiaofan/bootstrap probability-xiaofan/2.model"
datas <- dir(workdir)
datas <- datas[-grep(".txt|TPM|WPS-divide-bgWPS|WPS-divide-COV|WPS-substract-bgWPS|WPS-subtractWPS-divideCOV|miRNA",datas)]
comparisons <- c("CRC_HD","STAD_HD","STAD_CRC")
file_prefix <- "result_predict_FDR_top50"
nfiles <- 10
output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220629_xiaofan/bootstrap probability-xiaofan/summary"
dir.create(output_dir)

for( comparison in comparisons){
  #read.csv(paste0(workdir,"/",data,"/","CRC_HD"))
  example <- read.csv(paste0(workdir,"/",datas[1],"/",comparison,"/",file_prefix,"_1.txt"),sep = "\t")
  probability_all <- as.data.frame(matrix(numeric(0),ncol=0,
                                            nrow = nrow(example)))
  probability_all$ID <- example$ID
  probability_all$label <- example$label
  
  for(data in datas){
    
    probability_final <- as.data.frame(matrix(numeric(0),ncol=nfiles,
                                              nrow = nrow(example)))
    colnames(probability_final) <- c(1:nfiles)
    for(n in c(1:nfiles)){
      probability <- read.csv(paste0(workdir,"/",data,"/",comparison,"/",file_prefix,"_",n,".txt"),sep = "\t")
      probability_final[,n] <- probability[[paste0("test_",n)]]
    }
    probability_mean <- data.frame(data=rowMeans(probability_final))
    colnames(probability_mean) <- data
    probability_all <- cbind(probability_all,probability_mean)
    }
  write.table(probability_all,paste0(output_dir,"/",comparison,".txt"),sep = "\t", quote = FALSE,row.names = FALSE)
}

A <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220629_xiaofan/bootstrap probability-xiaofan/summary/STAD_CRC.txt",sep = "\t")
A <- as.matrix(A)
cvfit <- cv.glmnet(A[,-c(1,2)], A[,2] , standardize = TRUE , type.measure = "deviance", family = "binomial" , nfolds = 5, alpha = 1) 
cvfit
cvfit$lambda.min
cvfit$lambda.1se
coefficient <- coef(cvfit , s = "lambda.min")
#coefficient <- coef(cvfit, s = "lambda.1se")

Lasso_selected_genes <- as.data.frame(A[,(coefficient@i+1)[-1]])
genes_coefficients <- coefficient@x[-1]
