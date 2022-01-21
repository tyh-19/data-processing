#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(PRROC))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))

parser <- ArgumentParser(description='Merge model')
parser$add_argument('-i', '--input_dir', type='character', required=TRUE,
    help='input directory, which contains results like ./Expression/STADvsNC/ from LOO_v3_yuhuan.R.')
parser$add_argument('-f', '--AlterationNumber', type='integer', required=TRUE,
    help='Alteration number in final model.')
parser$add_argument('-p', '--positive', type='character', required=TRUE,
    choices=c('CRC','STAD','NC'))
parser$add_argument('-n','--negative', type='character', required=TRUE,
    choices=c('CRC','STAD','NC'))
parser$add_argument('-o', '--output_dir', type='character', required=TRUE,
    help='output file')
args <- parser$parse_args()

  input_dir <- args$input_dir
  output_dir_initial <- args$output_dir
  dir.create(output_dir_initial)
  positive <- args$positive
  negative <- args$negative
  n=args$AlterationNumber

  full_alterations <- dir(input_dir)
  #full_alterations <- full_alterations[full_alterations!="Merged"]
  full_alterations <- full_alterations[-grep("Merged",full_alterations)]
  full_alterations <- full_alterations[-grep("txt",full_alterations)]
  full_alterations <- full_alterations[full_alterations!="CNV"]
  full_alterations <- full_alterations[full_alterations!="Splicing_FDR0.05"]

  example <- read.csv(paste0(input_dir,"/",full_alterations[1],"/",positive,"vs",negative,"/",1,"_matrix_nolabel_selected.txt"),header = TRUE,sep = "\t")
  
  #cl <- makeCluster(16) #not to overload your computer
  #registerDoParallel(cl)

  k=1
  comb_alterations <- as.data.frame(combn(full_alterations,n))
  matrics_record <- as.data.frame(matrix(numeric(0),ncol=7,nrow=ncol(comb_alterations)))
  colnames(matrics_record) <- c("feature_number","combination","AUC","Specificity_0.95","Sensitivity_0.95","Specificity_youden","Sensitivity_youden")
  while(k<=ncol(comb_alterations)){
  alterations <- as.character(unlist(comb_alterations[k]))
  combination_name <- paste(paste(alterations,collapse = "|"),n,sep = "_")
  output_dir <- paste0(output_dir_initial,"/",combination_name)
  dir.create(output_dir)
  i=1
  prob_final <- {}
  sample_number <- nrow(example)
  while(i<=sample_number){ 
  j=1
  while(j<=length(alterations)){
    if(j==1){
    first <- read.csv(paste0(input_dir,"/",alterations[j],"/",positive,"vs",negative,"/",i,"_matrix_nolabel_selected.txt"),header = TRUE,sep = "\t")
    first$sample <- rownames(first)
    j=j+1
    } else { 
    second <- read.csv(paste0(input_dir,"/",alterations[j],"/",positive,"vs",negative,"/",i,"_matrix_nolabel_selected.txt"),header = TRUE,sep = "\t")
    second$sample <- rownames(second)
    first <- left_join(first,second,by=c("sample"="sample"))
    j=j+1
    }
  }
  
  rownames(first) <- first$sample
  first <- first[,-which(colnames(first)=="sample")]
  
  matrix <- first
  
  matrix$label <- unlist(lapply(strsplit(rownames(matrix),".",fixed = TRUE),function(x) x[1]))
  labels <- as.factor(matrix$label)
  
  matrix_nolabel <- matrix[,-which(colnames(matrix)=="label")]
  matrix_nolabel[is.na(matrix_nolabel)] <- 0
  matrix_nolabel_numeric <- mutate_all(matrix_nolabel, function(x) as.numeric(as.character(x)))
  rownames(matrix_nolabel_numeric) <- rownames(matrix_nolabel)
  
  #no feature selection
  matrix_nolabel_selected <- matrix_nolabel_numeric
  
  set.seed(6)
  train_labels <- as.factor(matrix[-i,]$label)
  #fix model parameters, such as mtry, ntree
  #tune_mtry <- tuneRF(matrix_nolabel_selected[-i,],train_labels) #tune mtry
  #tune_mtry <- as.data.frame(tune_mtry)
  #mtry_best <- tune_mtry[order(tune_mtry$OOBError,decreasing=FALSE),][1,1]
  
  #no mtry selection
  mtry_best <- floor(sqrt(ncol(matrix_nolabel_selected))) 
 
  #train RF model use n_sample-1
  RF_model <- randomForest(matrix_nolabel_selected[-i,],train_labels, mtry = mtry_best)
  
  #predict leaved 1 sample by trained model 
  test_labels <- factor(matrix[i,]$label,level=c(positive,negative))
  predict_prob <- predict(RF_model,newdata = matrix_nolabel_selected[i,], type = "prob")
  prob_temp <- as.data.frame(predict_prob)
  prob_final <- rbind(prob_final,prob_temp)
  i=i+1
  }

  write.table(prob_final,paste0(output_dir,"/prob_final_LOO.txt"),quote = FALSE,sep="\t")


  #report ROC
  predicted <- prob_final
  roc.curve <- roc(labels,predicted[,which(colnames(predicted)==positive)],levels=c(negative,positive),direction=c("<"))
  ci.auc(roc.curve,conf.level = 0.95)
  record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
  write.table(roc.curve$auc,paste0(output_dir,"/AUC.txt"),quote = FALSE,sep="\t")
  write.table(record$specificity[2],paste0(output_dir,"/Specificity_youden.txt"),quote = FALSE,sep="\t")
  write.table(record$recall[2],paste0(output_dir,"/Sensitivity_youden.txt"),quote = FALSE,sep="\t")

  #plot
  predicted$labels <- unlist(lapply(strsplit(rownames(predicted),".",fixed = TRUE),function(x) x[1]))
  predicted$labels <- gsub(positive,1,predicted$labels)
  predicted$labels <- gsub(negative,0,predicted$labels)
  predicted$labels <- as.numeric(predicted$labels)

  pdf(paste0(output_dir,"/AUROC.pdf"),width=7.8,height=6)
  roc_curve_prediction_df <- PRROC::roc.curve(scores.class0=predicted[,which(colnames(predicted)==positive)], weights.class0 = predicted$labels,
                                            curve=TRUE,rand.compute = TRUE) # rand.compute 是0.5的随机分类器，为了画对角线的

  par(cex.axis=2, cex.lab =1.6,font.axis =1.6, mar = c(5,5,2,3))
  plot(roc_curve_prediction_df,rand.plot = TRUE,auc.main=FALSE,legend =4,lwd=5)
  text(1.2, 0.4, paste('AUROC =', ' ',round(roc_curve_prediction_df$auc, digits = 4), sep =""), cex=2, font=2, col='black')
  dev.off()

  pdf(paste0(output_dir,"/AUPR.pdf"),width=7.8,height=6)
  pr_curve_prediction_df <- PRROC::pr.curve(scores.class0=predicted[,which(colnames(predicted)==positive)], weights.class0 = predicted$labels,
                                          curve=TRUE) # rand.compute 是0.5的随机分类器，为了画对角线的

  par(cex.axis=2, cex.lab =1.6,font.axis =1.6, mar = c(5,5,2,3))
  plot(pr_curve_prediction_df,rand.plot = TRUE,auc.main=FALSE,legend =4,lwd=5)
  text(1.2, 0.4, paste('AUPR =', ' ',round(pr_curve_prediction_df$auc.integral, digits = 4), sep =""), cex=2, font=2, col='black')
  dev.off()

  pred <- predicted[,which(colnames(predicted)==positive)]
  ref <- unlist(lapply(strsplit(rownames(predicted),".",fixed = TRUE),function(x) x[1]))

  m=1
  #specificity 0.95 threhold
  threhold <- sort(predicted[which(predicted$labels=="0"),which(colnames(predicted)==positive)])[ceiling(0.95*nrow(predicted[which(predicted$labels=="0"),]))]
  while(m<=length(pred)){
  if(pred[m]>threhold){
    pred[m] <- positive
  } else {
    pred[m] <- negative
  }
    m=m+1
  }

  pred <- factor(pred,levels = c(positive,negative))
  ref <- factor(ref,levels = c(positive,negative))

  confusion_matrix <- confusionMatrix(pred,ref)
  write.table(confusion_matrix$table,paste0(output_dir,"/confusion_matrix_0.95specificity.txt"),quote = FALSE,sep="\t")
  write.table(confusion_matrix$byClass[1],paste0(output_dir,"/Sensitivity_0.95specificity.txt"),quote = FALSE,sep="\t")
  write.table(confusion_matrix$byClass[2],paste0(output_dir,"/Specificity_0.95specificity.txt"),quote = FALSE,sep="\t")
  write.table(confusion_matrix$byClass[7],paste0(output_dir,"/F1_0.95specificity.txt"),quote = FALSE,sep="\t")
  write.table(confusion_matrix$byClass[11],paste0(output_dir,"/Balanced_accuracy_0.95specificity.txt"),quote = FALSE,sep="\t")

  matrics_record[k,"feature_number"] <- n
  matrics_record[k,"combination"] <- combination_name
  matrics_record[k,"AUC"] <- roc.curve$auc
  matrics_record[k,"Specificity_0.95"] <- confusion_matrix$byClass[2]
  matrics_record[k,"Sensitivity_0.95"] <- confusion_matrix$byClass[1]
  matrics_record[k,"Specificity_youden"] <- record$specificity[2]
  matrics_record[k,"Sensitivity_youden"] <-record$recall[2]
  k=k+1  
  }
  write.table(matrics_record,paste0(output_dir_initial,"/",n,"_matrics_record.txt"),quote = FALSE,sep="\t")
  #stopCluster(cl)
