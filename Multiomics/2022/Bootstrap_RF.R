
#! /usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(PRROC))
suppressPackageStartupMessages(library(dplyr))

parser <- ArgumentParser(description='Bootstrap Random Forest Model')
parser$add_argument('-s', '--samples', type='character', required=TRUE,
                    help='input samples and labels. Column names: sample_id, Group, miRNA_id, RNA_id, DNA_id, Methylation_id.')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input count matrix. Rows are features. Columns are samples.')
parser$add_argument('-t', '--sample_type', type='character', required=TRUE,
                    choices=c("DNA_id","Methylation_id","RNA_id","miRNA_id"))
parser$add_argument('-r', '--repeat_number', type='integer', required=TRUE,
                    help='repeat number')
parser$add_argument('-bsn', '--bootstrap_sample_number', type='integer', required=TRUE,
                    help='bootstrap sample number')
parser$add_argument('-b', '--bootstrap_number', type='integer', required=TRUE,
                    help='bootstrap times')
parser$add_argument('-f', '--featureNumber', type='integer', required=TRUE,
                    help='feature number in final model. 200 is suggested. If total features are less than 200, use all features.')
parser$add_argument('-p', '--positive', type='character', required=TRUE,
                    choices=c('CRC','STAD','HD'))
parser$add_argument('-n','--negative', type='character', required=TRUE,
                    choices=c('CRC','STAD','HD'))
parser$add_argument('-o', '--outdir', type='character', required=TRUE,
                    help='output directory (do not end with /)')
args <- parser$parse_args()


message('read count matrix: ', args$matrix)
samples <- read.csv(args$sample,row.names = NULL)
data <- read.table(args$matrix,sep = "\t",check.names = FALSE,header = TRUE, row.names = 1)
sample_type <- args$sample_type #c("DNA_id","Methylation_id","RNA_id","miRNA_id")
positive_label <- args$positive
negative_label <- args$negative
bootstrap_sample_number <- args$bootstrap_sample_number
feature_number <- args$featureNumber
bootstrap_number <- args$bootstrap_number
repeat_number <- args$repeat_number
outdir <- paste0(args$outdir,"/",positive_label,"_vs_",negative_label)

#Bootstrap random forest
#preparation
library(randomForest)
library(caret)
#function
#commen function for differential analysis: wilcox test
Wilcox_test <- function(mat_raw,des,output_res){
  norm_method <- 'NA'
  samples <- des$samples
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
  
  #filter events < 20% samples in minial group
  #positive_prop <- rowSums(matrix[,positive] > 0)/length(positive)
  #negative_prop <- rowSums(matrix[,negative] > 0)/length(negative)
  #matrix <- matrix[(negative_prop>=0.2)&(positive_prop>=0.2),]
  
  #write.table(matrix, output_matrix, sep='\t', quote=FALSE, row.names=TRUE)
  
  pvalues <- apply(matrix, 1, test_func)
  matrix_logcpm = log2(matrix + 1)
  logFC <- apply(matrix_logcpm[,positive], 1, mean) -
    apply(matrix_logcpm[,negative], 1, mean)
  deltaFraction <- apply(matrix[,positive], 1, mean) -
    apply(matrix[,negative], 1, mean)
  
  positive_gini <- as.numeric(gini(t(matrix[,positive])))
  negative_gini <- as.numeric(gini(t(matrix[,negative])))
  
  res <- data.frame(log2FoldChange=logFC,
                    deltaFraction=deltaFraction,
                    pvalue=pvalues, 
                    padj=p.adjust(pvalues, method='BH'),
                    baseMean=apply(matrix, 1, mean),
                    positive_gini = positive_gini,
                    negative_gini = negative_gini)
  res <- res[order(res$pvalue,decreasing = FALSE),]
  write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE)
  res
}

#input: sample, data, sample_type, negative_label, positive label
samples <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/ML_samples.csv",row.names = NULL)
data <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Alternative promoter/matrix_forML/Altpromoter_ML.txt",sep = "\t",check.names = FALSE,header = TRUE, row.names = 1)
sample_type <- "RNA_id" #c("DNA_id","Methylation_id","RNA_id","miRNA_id")
positive_label <- "CRC"
negative_label <- "HD"
bootstrap_sample_number <- 50#50
feature_number <- 3000#3000
bootstrap_number <- 50#5
repeat_number <- 3#3
outdir <- paste0("/Users/yuhuan/Desktop/Alternative promoter/",positive_label,"_vs_",negative_label)
dir.create(outdir, recursive = TRUE)

samples <- samples[complete.cases(samples[[sample_type]]),]

positive_samples <- samples[grep(positive_label,samples$Group),]

negative_samples <- samples[grep(negative_label,samples$Group),]

ML_samples <- rbind(positive_samples,negative_samples)
#bootstrap_sample_number <- nrow(ML_samples)-1

j=1
aggregated_auc_recall_specificity <- as.data.frame(matrix(numeric(0),ncol=4))
colnames(aggregated_auc_recall_specificity) <- c("Repeat","Aggregated_AUC","Aggregated_recall","Aggregated_specificity")

while(j <= repeat_number){
i=1
prob_train_final <- {}
prob_final <- {}
test_auc_recall_specificity <- as.data.frame(matrix(numeric(0),ncol=4))
colnames(test_auc_recall_specificity) <- c("Random_seed","AUC","Recall","Specificity")
while(i <= bootstrap_number){
random_seed <- 100000+1000*j+i
set.seed(random_seed)
index <- sample(nrow(ML_samples), size = bootstrap_sample_number, replace = TRUE)
train_samples <- as.data.frame(ML_samples[index,])
test_samples <- as.data.frame(ML_samples[-unique(index),])

#data split
train_data <- as.data.frame(t(data[,as.character(train_samples[[sample_type]])]))
train_data_forDE <- as.data.frame(data[,as.character(train_samples[[sample_type]])])
#train_data[is.na(train_data)] <- 0
#train_data_forDE[is.na(train_data_forDE)] <- 0
train_label <- as.character(train_samples$Group)
train_label <- factor(train_label, levels = c(positive_label,negative_label))

test_data <- as.data.frame(t(data[,as.character(test_samples[[sample_type]])]))
colnames(test_data) <- rownames(data)
rownames(test_data) <- as.character(test_samples[[sample_type]])
#test_data[is.na(test_data)] <- 0
test_label <- as.character(test_samples$Group)
test_label <- as.character(test_samples$Group)
test_label <- factor(test_label, levels = c(positive_label,negative_label))

#feature selection
des <- data.frame("samples"=rownames(train_data),"group"=train_label)
des$samples <- as.character(des$samples)
des$group <- gsub(positive_label,"positive",des$group)
des$group <- gsub(negative_label,"negative",des$group)
output_res <- paste0(outdir,"/feature_selected/",j,"_",random_seed,".txt")
dir.create(paste0(outdir,"/feature_selected/"))
mat_raw <- train_data_forDE
Differential_result <- Wilcox_test(mat_raw,des,output_res)
selected_features <- rownames(head(Differential_result,feature_number))
train_data <- train_data[,selected_features]
train_data[is.na(train_data)] <- 0

#model training (by caret)
#tune_mtry <- tuneRF(train_data,train_label) #tune mtry
need_tune <- TRUE
tune_times=1
while(need_tune){
repeat_tune <- FALSE
tryCatch(tune_mtry <- tuneRF(train_data,train_label), error = function(e) {
  repeat_tune <<- TRUE
  message("Failed in first tuneRF, try again.")
  })
if(repeat_tune) {need_tune <- TRUE
tune_times=tune_times+1} else {need_tune <- FALSE}
}
message("Tune parameter times: ",tune_times)

#if(repeat_tune) { 
#  message("Failed in first tuneRF, try again.")
#  tune_mtry <- tuneRF(train_data,train_label) 
#  }
tune_mtry <- as.data.frame(tune_mtry)
mtry_best <- tune_mtry[order(tune_mtry$OOBError,decreasing=FALSE),][1,1]

#train RF model use n_sample-1? or setted bootstrap_sample_number
RF_model <- randomForest(train_data,train_label, mtry = mtry_best)

train_predicted_tmp <- as.data.frame(RF_model$votes)
train_roc.curve <- roc(train_label,train_predicted_tmp[,which(colnames(train_predicted_tmp)==positive_label)],levels=c(negative_label,positive_label),direction=c("<"))

train_predicted_tmp$seed <- random_seed
train_predicted_tmp$samples <- rownames(train_predicted_tmp)
train_predicted_tmp$samples <- as.character(lapply(strsplit(x=train_predicted_tmp$samples,split = ".",fixed = TRUE), function(x) x[1]))
train_predicted_tmp$label <- train_label
rownames(train_predicted_tmp) <- NULL
prob_train_final <- rbind(prob_train_final,train_predicted_tmp)

#testing
#predict unsampled samples by trained model
test_data <- test_data[,selected_features]
test_data[is.na(test_data)] <- 0
predict_prob <- predict(RF_model,newdata = test_data, type = "prob")
test_roc.curve <- roc(test_label,predict_prob[,which(colnames(predict_prob)==positive_label)],levels=c(negative_label,positive_label),direction=c("<"))
test_AUC <- test_roc.curve$auc[1]

test_record <- ci.coords(test_roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
test_recall <- test_record$recall[2]
test_specificity <- test_record$specificity[2]
test_auc_recall_specificity[i,1] <- random_seed
test_auc_recall_specificity[i,2] <- test_AUC
test_auc_recall_specificity[i,3] <- test_recall
test_auc_recall_specificity[i,4] <- test_specificity 


#test_label_numeric <- gsub(positive_label,"1",test_label)
#test_label_numeric <- gsub(negative_label,"0",test_label_numeric)
#test_label_numeric <- as.numeric(test_label_numeric)
#roc_curve_prediction_df <- PRROC::roc.curve(scores.class0=predict_prob[,which(colnames(predict_prob)==positive_label)], weights.class0 = test_label_numeric,
#                                            curve=TRUE,rand.compute = TRUE) # rand.compute 是0.5的随机分类器，为了画对角线的

#par(cex.axis=2, cex.lab =1.6,font.axis =1.6, mar = c(5,5,2,3))
#plot(roc_curve_prediction_df,rand.plot = TRUE,auc.main=FALSE,legend =4,lwd=5)
#text(1.3, 0.2, paste('AUROC =', ' ',round(roc_curve_prediction_df$auc, digits = 4), sep =""), cex=2, font=2, col='black')
#dev.off()

prob_temp <- as.data.frame(predict_prob)
prob_temp$seed <- random_seed
prob_temp$samples <- rownames(prob_temp)
prob_temp$label <- test_label
rownames(prob_temp) <- NULL
prob_final <- rbind(prob_final,prob_temp)

i=i+1
}

#all train prob
output_train_prob <- paste0(outdir,"/train_probs/",j,"_train_probability_of_each_sampling.txt")
dir.create(paste0(outdir,"/train_probs"))
write.table(prob_train_final,output_train_prob,quote = FALSE,sep="\t",row.names = FALSE)

#all test prob
output_test_prob <- paste0(outdir,"/test_probs/",j,"_test_probability_of_each_sampling.txt")
dir.create(paste0(outdir,"/test_probs"))
write.table(prob_final,output_test_prob,quote = FALSE,sep="\t",row.names = FALSE)

#all test performance
output_test_performance <- paste0(outdir,"/test_performance/",j,"_AUC.txt")
dir.create(paste0(outdir,"/test_performance"))
write.table(test_auc_recall_specificity,output_test_performance,quote = FALSE,sep="\t",row.names = FALSE)

#aggregated performance
prob_final$label <- gsub(positive_label,"1",prob_final$label)
prob_final$label <- gsub(negative_label,"0",prob_final$label)
prob_final$label <- as.numeric(prob_final$label)
prob_final_aggregated <- aggregate(prob_final,by = list("samples"=prob_final$samples),FUN = mean)
prob_final_aggregated$label <- gsub("1",positive_label,prob_final_aggregated$label)
prob_final_aggregated$label <- gsub("0",negative_label,prob_final_aggregated$label)

aggregated_AUC <- roc(prob_final_aggregated$label,prob_final_aggregated[,which(colnames(prob_final_aggregated)==positive_label)],levels=c(negative_label,positive_label),direction=c("<"))
aggregated_record <- ci.coords(aggregated_AUC,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
aggregated_recall <- aggregated_record$recall[2]
aggregated_specificity <- aggregated_record$specificity[2]

aggregated_auc_recall_specificity[j,1] <- j
aggregated_auc_recall_specificity[j,2] <- aggregated_AUC$auc[1]
aggregated_auc_recall_specificity[j,3] <- aggregated_recall
aggregated_auc_recall_specificity[j,4] <- aggregated_specificity


output_aggregated_performance <- paste0(outdir,"/aggregated_performance/aggregated_AUC.txt")
dir.create(paste0(outdir,"/aggregated_performance"))
write.table(aggregated_auc_recall_specificity,output_aggregated_performance,quote = FALSE,sep="\t",row.names = FALSE)

j=j+1
}



