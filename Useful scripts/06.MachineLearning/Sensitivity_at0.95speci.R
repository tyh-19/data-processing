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
parser$add_argument('-i', '--input_dir', type='character', required=TRUE,
    help='input dir contains a dataframe named prob_final_LOO.txt. colnames are positive and negative group. rownames are samples which named "positive.xxx"/"negative.xxx".')
parser$add_argument('-p','--positive', type='character', required=TRUE,
    choices=c('CRC','STAD','NC'))
parser$add_argument('-n','--negative', type='character', required=TRUE,
    choices=c('CRC','STAD','NC'))
parser$add_argument('-o', '--output_dir', type='character', required=TRUE,
    help='output file')
args <- parser$parse_args()

input_dir <- args$input_dir
output_dir <- args$output_dir
positive <- args$positive
negative <- args$negative

prob_final <- read.csv(paste0(input_dir,"/prob_final_LOO.txt"),header = TRUE, row.names = 1, sep="\t")

#report ROC
predicted <- prob_final
ref <- unlist(lapply(strsplit(rownames(predicted),".",fixed = TRUE),function(x) x[1]))
roc.curve <- roc(ref,predicted[,which(colnames(predicted)==positive)],levels=c(negative,positive),direction=c("<"))
ci.auc(roc.curve,conf.level = 0.95)
record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
write.table(roc.curve$auc,paste0(output_dir,"/AUC_recalculated.txt"),quote = FALSE,sep="\t")
write.table(record$specificity[2],paste0(output_dir,"/Specificity_youden.txt"),quote = FALSE,sep="\t")
write.table(record$recall[2],paste0(output_dir,"/Sensitivity_youden.txt"),quote = FALSE,sep="\t")

#plot
predicted$labels <- unlist(lapply(strsplit(rownames(predicted),".",fixed = TRUE),function(x) x[1]))
predicted$labels <- gsub(positive,1,predicted$labels)
predicted$labels <- gsub(negative,0,predicted$labels)
predicted$labels <- as.numeric(predicted$labels)

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

