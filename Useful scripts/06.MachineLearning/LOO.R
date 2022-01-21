#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(PRROC))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
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

#remove highly correlated columne
matrix_cor = cor(matrix_nolabel_numeric)
highly_correlated = findCorrelation(matrix_cor, cutoff=0.9, exact = TRUE)
highly_correlated = sort(highly_correlated)
matrix_nolabel_numeric = matrix_nolabel_numeric[,-c(highly_correlated)]
write.table(matrix_nolabel_numeric,paste0(outdir,"matrix_nolabel_rm_higly_correlated.txt"),quote = FALSE,sep="\t")
write.table(higly_correlated,paste0(outdir,"highly_correlated_features.txt"),quote = FALSE,sep="\t")

labels <- as.factor(matrix$label)

#feature selection top200
rfFuncs$summary <- twoClassSummary

set.seed(123)
subsetSizes <- c(featureNumber)
seeds <- vector(mode = "list", length = nrow(matrix_nolabel_numeric)+1 )
for(i in 1:(nrow(matrix_nolabel_numeric))) seeds[[i]] <- sample.int(1000, length(subsetSizes) + 1)
seeds[[(nrow(matrix_nolabel_numeric)+1)]] <- sample.int(1000, 1)

set.seed(1)
rfectrl <- rfeControl(functions=rfFuncs,
                      seeds = seeds,
                      saveDetails = TRUE,
                      verbose = TRUE,
                      method="LOOCV")

rfe.results <- rfe(matrix_nolabel_numeric,labels,
                   sizes = c(featureNumber),
                   rfeControl = rfectrl,
                   metric = "ROC")

#get feature importance
y <- rfe.results$variables
finalImp <- ddply(y[, c("Overall", "var")], .(var), function(x) mean(x$Overall,na.rm = TRUE))
names(finalImp)[2] <- "Overall"
finalImp <- finalImp[order(finalImp$Overall, decreasing = TRUE),]

if(rfe.results$optsize==featureNumber) {
  message(paste0(featureNumber," is RFE suggested. Using top"))
  feature_selected <- predictors(rfe.results)
} else if(featureNumber<=ncol(matrix_nolabel_numeric)) {
  message(paste0(featureNumber," is not suggested. Using top",featureNumber," as final feature."))
  feature_selected <- as.character(finalImp$var[1:featureNumber])
} else {
  total_feature_number <- ncol(matrix_nolabel_numeric)
  message(paste0(featureNumber," is larger than total feature number. Using total feature number:",total_feature_number," as final feature."))
  feature_selected <- as.character(finalImp$var[1:total_feature_number])
}

matrix_nolabel_selected <- matrix_nolabel_numeric[,feature_selected]
write.table(matrix_nolabel_selected,paste0(outdir,"matrix_nolabel_selected.txt"),quote = FALSE,sep="\t")
write.table(finalImp,paste0(outdir,"feature_importance.txt"),quote = FALSE,sep="\t")

i=1
prob_final <- {}
while(i<=nrow(matrix)) {
train_labels <- as.factor(matrix[-i,]$label)
#fix model parameters, such as mtry, ntree
tune_mtry <- tuneRF(matrix_nolabel_selected[-i,],train_labels) #tune mtry
tune_mtry <- as.data.frame(tune_mtry)
mtry_best <- tune_mtry[order(tune_mtry$OOBError,decreasing=FALSE),][1,1]

#train RF model use n_sample-1
RF_model <- randomForest(matrix_nolabel_selected[-i,],train_labels, mtry = mtry_best)

#predict leaved 1 sample by trained model 
test_labels <- factor(matrix[i,]$label,level=c(positive,negative))
predict_prob <- predict(RF_model,newdata = matrix_nolabel_selected[i,], type = "prob")
prob_temp <- as.data.frame(predict_prob)
prob_final <- rbind(prob_final,prob_temp)
i=i+1
}
write.table(prob_final,paste0(outdir,"prob_final_LOO.txt"),quote = FALSE,sep="\t")


#report ROC
predicted <- prob_final
roc.curve <- roc(labels,predicted[,which(colnames(predicted)==positive)])
ci.auc(roc.curve,conf.level = 0.95)
record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
write.table(record,paste0(outdir,"matrices_record.txt"),quote = FALSE,sep="\t")

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

confusionMatrix <- confusionMatrix(pred,ref)
write.table(confusionMatrix,paste0(outdir,"confusion_matrix.txt"),quote = FALSE,sep="\t")

