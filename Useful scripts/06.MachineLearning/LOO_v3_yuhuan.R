#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(PRROC))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(e1071))
parser <- ArgumentParser(description='Random Forest LOO')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('-f', '--featureNumber', type='integer', required=TRUE,
    help='feature number in final model. 200 is suggested. If total features are less than 200, use all features.')
parser$add_argument('-p', '--positive', type='character', required=TRUE,
    choices=c('CRC','STAD','NC'))
parser$add_argument('-n','--negative', type='character', required=TRUE,
    choices=c('CRC','STAD','NC'))
parser$add_argument('-o', '--outdir', type='character', required=TRUE,
    help='output file')
args <- parser$parse_args()

  cl <- makeCluster(16) #not to overload your computer
  registerDoParallel(cl)

message('read count matrix: ', args$matrix)
outdir <- args$outdir
positive <- args$positive
negative <- args$negative
featureNumber <- args$featureNumber
#trasnform data matrix, col are features and row are samples, then add class/label at last column
raw_matrix <- read.csv(args$matrix,sep = "\t",header = TRUE, row.names = 1)
full_matrix <- as.data.frame(t(raw_matrix))
full_matrix$label <- unlist(lapply(strsplit(rownames(full_matrix),".",fixed = TRUE),function(x) x[1]))
  
  positive_matrix <- full_matrix[which(full_matrix$label==positive),]
  negative_matrix <- full_matrix[which(full_matrix$label==negative),]
  matrix <- rbind(positive_matrix,negative_matrix)
  #numeric matrix
  matrix_nolabel <- matrix[,-which(colnames(matrix)=="label")]
  matrix_nolabel[is.na(matrix_nolabel)] <- 0
  matrix_nolabel_numeric <- mutate_all(matrix_nolabel, function(x) as.numeric(as.character(x)))
  rownames(matrix_nolabel_numeric) <- rownames(matrix_nolabel)
  
  labels <- as.factor(matrix$label)
  
  i=1
  prob_final <- {}
  train_auc <- as.data.frame(matrix(numeric(0),ncol=2))
  while(i<=nrow(matrix)) {
    train_labels <- as.factor(matrix[-i,]$label)
    
    #select features in training set
    {
      #feature selection top200
      rfFuncs$summary <- twoClassSummary
      
      #set.seed(-100000-i)
      subsetSizes <- c(featureNumber)
      seeds <- vector(mode = "list", length = nrow(matrix_nolabel_numeric[-i,])+1 )
      k=1
      #for(k in 1:(nrow(matrix_nolabel_numeric[-k,]))) seeds[[k]] <- sample.int(100000, length(subsetSizes) + 1)
      #seeds[[(nrow(matrix_nolabel_numeric[-k,])+1)]] <- sample.int(100000, 1)
      for(k in 1:(nrow(matrix_nolabel_numeric[-k,]))) seeds[[k]] <- c(i, i + 1)
      seeds[[(nrow(matrix_nolabel_numeric[-k,])+1)]] <- sample.int(i+2, 1)
      
      set.seed(-i)
      rfectrl <- rfeControl(functions=rfFuncs,
                            seeds = seeds,
                            saveDetails = TRUE,
                            verbose = TRUE,
                            method="LOOCV")
      
      rfe.results <- rfe(matrix_nolabel_numeric[-i,],train_labels,
                         sizes = c(featureNumber),
                         rfeControl = rfectrl,
                         metric = "ROC")
      
      #get feature importance
      y <- rfe.results$variables
      finalImp <- ddply(y[, c("Overall", "var")], .(var), function(x) mean(x$Overall,na.rm = TRUE))
      names(finalImp)[2] <- "Overall"
      finalImp <- finalImp[order(finalImp$Overall, decreasing = TRUE),]
      
      if(rfe.results$optsize==featureNumber) {
        print(paste0(featureNumber," is RFE suggested. Using top"))
        feature_selected <- predictors(rfe.results)
      } else if(featureNumber<=ncol(matrix_nolabel_numeric[-i,])){
        print(paste0(featureNumber," is not suggested. Using top",featureNumber," as final feature."))
        feature_selected <- as.character(finalImp$var[1:featureNumber])
      } else {
        total_feature_number <- ncol(matrix_nolabel_numeric[-i,])
        print(paste0(featureNumber," is larger than total feature number. Using total feature number:",total_feature_number," as final feature."))
        feature_selected <- as.character(finalImp$var[1:total_feature_number])
      }
      
      
      matrix_nolabel_selected <- matrix_nolabel_numeric[,feature_selected]
      write.table(matrix_nolabel_selected,paste0(outdir,i,"_","matrix_nolabel_selected.txt"),quote = FALSE,sep="\t")
      write.table(finalImp,paste0(outdir,i,"_","feature_importance.txt"),quote = FALSE,sep="\t")
    }
    
    set.seed(100000+i)
    #fix model parameters, such as mtry, ntree
    tune_mtry <- tuneRF(matrix_nolabel_numeric[-i,],train_labels) #tune mtry
    tune_mtry <- as.data.frame(tune_mtry)
    mtry_best <- tune_mtry[order(tune_mtry$OOBError,decreasing=FALSE),][1,1]
    
    #train RF model use n_sample-1
    RF_model <- randomForest(matrix_nolabel_selected[-i,],train_labels, mtry = mtry_best)
    
    train_predicted <- as.data.frame(RF_model$votes)
    train_roc.curve <- roc(train_labels,train_predicted[,which(colnames(train_predicted)==positive)],levels=c(negative,positive),direction=c("<"))
    
    #train_auc
    train_auc[i,1] <- i
    train_auc[i,2] <- train_roc.curve$auc
    write.table(train_auc,paste0(outdir,"train_auc.txt"),quote = FALSE,sep="\t")
    
    #predict leaved 1 sample by trained model 
    test_labels <- factor(matrix[i,]$label,level=c(positive,negative))
    predict_prob <- predict(RF_model,newdata = matrix_nolabel_selected[i,], type = "prob")
    prob_temp <- as.data.frame(predict_prob)
    prob_final <- rbind(prob_final,prob_temp)
    i=i+1
  }
  write.table(train_auc,paste0(outdir,"train_auc.txt"),quote = FALSE,sep="\t")
  write.table(prob_final,paste0(outdir,"prob_final_LOO.txt"),quote = FALSE,sep="\t")
  
  
  #report ROC
  predicted <- prob_final
  roc.curve <- roc(labels,predicted[,which(colnames(predicted)==positive)],levels=c(negative,positive),direction=c("<"))
  ci.auc(roc.curve,conf.level = 0.95)
  record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
  write.table(roc.curve$auc,paste0(outdir,"AUC.txt"),quote = FALSE,sep="\t")
  #write.table(record,paste0(outdir,"matrices_record.txt"),quote = FALSE,sep="\t")
  
  #plot
  predicted$labels <- unlist(lapply(strsplit(rownames(predicted),".",fixed = TRUE),function(x) x[1]))
  predicted$labels <- gsub(positive,1,predicted$labels)
  predicted$labels <- gsub(negative,0,predicted$labels)
  predicted$labels <- as.numeric(predicted$labels)
  
  pdf(paste0(outdir,"AUROC.pdf"),width=7.8,height=6)
  roc_curve_prediction_df <- PRROC::roc.curve(scores.class0=predicted[,which(colnames(predicted)==positive)], weights.class0 = predicted$labels,
                                              curve=TRUE,rand.compute = TRUE) # rand.compute 是0.5的随机分类器，为了画对角线的
  
  par(cex.axis=2, cex.lab =1.6,font.axis =1.6, mar = c(5,5,2,3))
  plot(roc_curve_prediction_df,rand.plot = TRUE,auc.main=FALSE,legend =4,lwd=5)
  text(1.3, 0.2, paste('AUROC =', ' ',round(roc_curve_prediction_df$auc, digits = 4), sep =""), cex=2, font=2, col='black')
  dev.off()
  
  pdf(paste0(outdir,"AUPR.pdf"),width=7.8,height=6)
  pr_curve_prediction_df <- PRROC::pr.curve(scores.class0=predicted[,which(colnames(predicted)==positive)], weights.class0 = predicted$labels,
                                            curve=TRUE) # rand.compute 是0.5的随机分类器，为了画对角线的
  
  par(cex.axis=2, cex.lab =1.6,font.axis =1.6, mar = c(5,5,2,3))
  plot(pr_curve_prediction_df,rand.plot = TRUE,auc.main=FALSE,legend =4,lwd=5)
  text(1.3, 0.2, paste('AUPR =', ' ',round(pr_curve_prediction_df$auc.integral, digits = 4), sep =""), cex=2, font=2, col='black')
  dev.off()
  
  pred <- predicted[,which(colnames(predicted)==positive)]
  ref <- unlist(lapply(strsplit(rownames(predicted),".",fixed = TRUE),function(x) x[1]))
  
  i=1
  while(i<=length(pred)){
    if(pred[i]>=0.5){
      pred[i] <- positive
    } else {
      pred[i] <- negative
    }
    i=i+1
  }
  
  pred <- factor(pred,levels = c(positive,negative))
  ref <- factor(ref,levels = c(positive,negative))
  
  confusion_matrix <- confusionMatrix(pred,ref)
  write.table(confusion_matrix$table,paste0(outdir,"confusion_matrix.txt"),quote = FALSE,sep="\t")
  write.table(confusion_matrix$byClass[1],paste0(outdir,"Sensitivity.txt"),quote = FALSE,sep="\t")
  write.table(confusion_matrix$byClass[2],paste0(outdir,"Specificity.txt"),quote = FALSE,sep="\t")
  write.table(confusion_matrix$byClass[7],paste0(outdir,"F1.txt"),quote = FALSE,sep="\t")
  write.table(confusion_matrix$byClass[11],paste0(outdir,"Balanced_accuracy.txt"),quote = FALSE,sep="\t")
  
  stopCluster(cl)
