library(glmnet)
library(pROC)
#Late integration of different modals/omics
workdir <- "/Users/yuhuan/Desktop/result_3000f_3b_50bsn/"
setwd(workdir)
datas <- dir(workdir)
repeat_times <- 3
comparisons <- c("CRC_vs_HD","STAD_vs_HD","STAD_vs_CRC")
outdir <- "/Users/yuhuan/Desktop/summary_test"
dir.create(outdir)

j=1
while(j<=repeat_times){
for(comparison in comparisons){
positive <- as.character(unlist(strsplit(comparison,"_vs_",fixed = TRUE))[1])
negative <- as.character(unlist(strsplit(comparison,"_vs_",fixed = TRUE))[2])
dir.create(paste0(outdir,"/",comparison))
train_probs_1 <- read.table(paste0(datas[1],"/",comparison,"/train_probs/",j,"_train_probability_of_each_sampling.txt"),sep = "\t",header = TRUE)
all_train_probs <- train_probs_1[,c("seed","samples","label")]
test_probs_1 <- read.table(paste0(datas[1],"/",comparison,"/test_probs/",j,"_test_probability_of_each_sampling.txt"),sep = "\t",header = TRUE)
all_test_probs <- test_probs_1[,c("seed","samples","label")]
test_performance_1 <- read.table(paste0(datas[1],"/",comparison,"/test_performance/",j,"_AUC.txt"),sep = "\t",header = TRUE)
all_test_performance <- test_performance_1[,c("Random_seed"),drop = FALSE]
all_test_performance_AUC <- all_test_performance
all_test_performance_spec_youden <- all_test_performance
all_test_performance_recall_youden <- all_test_performance
all_test_performance_precision_youden <- all_test_performance
all_test_performance_recall_100spec <- all_test_performance

for(data in datas){
train_probs <- read.table(paste0(data,"/",comparison,"/train_probs/",j,"_train_probability_of_each_sampling.txt"),sep = "\t",header = TRUE)
train_probs <- train_probs[,1,drop = FALSE]
colnames(train_probs) <- paste0(colnames(train_probs),"-",strsplit(data,".txt")[[1]][1])

test_probs <- read.table(paste0(data,"/",comparison,"/test_probs/",j,"_test_probability_of_each_sampling.txt"),sep = "\t",header = TRUE)
test_probs_forothermetric <- test_probs
test_probs <- test_probs[,1,drop = FALSE]
colnames(test_probs) <- paste0(colnames(test_probs),"-",strsplit(data,".txt")[[1]][1])

test_performance <- as.data.frame(matrix(numeric(0),ncol=7))
colnames(test_performance) <- c("comparison","seed","AUC","spec_youden","recall_youden","precision_youden","recall_100%spec")
for(seed in unique(test_probs_forothermetric$seed)){
  single_alteration_prob <- test_probs_forothermetric[test_probs_forothermetric$seed==seed,positive]
  single_alteration_label <- test_probs_forothermetric[test_probs_forothermetric$seed==seed,"label"]
  single_alteration_label <- gsub(positive,"1",single_alteration_label)
  single_alteration_label <- gsub(negative,"0",single_alteration_label)
  single_alteration_label <- as.numeric(single_alteration_label)
  #AUC
  test_roc.curve <- roc(single_alteration_label,single_alteration_prob,levels=c(0,1),direction=c("<"))
  test_AUC <- test_roc.curve$auc[1]
  print(test_AUC)
  #other metric
  record <- ci.coords(test_roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
  spec_youden <- record$specificity[2]
  recall_youden <-record$recall[2]
  precision_youden <-record$precision[2]
  #recall at 100% specificity
  ref <- gsub("1","positive",single_alteration_label)
  ref <- gsub("0","negative",ref)
  cutoff = 0.5
  cutoff = max(single_alteration_prob[ref=="negative"])
  i=1
  pred={}
  while(i<=length(single_alteration_prob)){
    if(single_alteration_prob[i]>cutoff){
      pred[i] <- "positive"
    } else {
      pred[i] <- "negative"
    }
    i=i+1
    
  }
  
  pred <- factor(pred,levels = c("positive","negative"))
  ref <- factor(ref,levels = c("positive","negative"))
  confusion_matrix <- confusionMatrix(pred,ref)
  
  recall_100 <- confusion_matrix$byClass[1]
  
  test_performance_tmp <- data.frame("comparison"=comparison,"seed"=seed,"AUC"=test_AUC,"spec_youden"=spec_youden,"recall_youden"=recall_youden,"precision_youden"=precision_youden,"recall_100%spec"=recall_100)
  test_performance <- rbind(test_performance,test_performance_tmp)
}

test_performance_AUC <- read.table(paste0(data,"/",comparison,"/test_performance/",j,"_AUC.txt"),sep = "\t",header = TRUE)
test_performance_AUC <- test_performance_AUC[,c("AUC"),drop = FALSE]
colnames(test_performance_AUC) <- paste0(colnames(test_performance_AUC),"-",strsplit(data,".txt")[[1]][1])

test_performance_spec_youden <- test_performance[,c("spec_youden"),drop = FALSE]
colnames(test_performance_spec_youden) <- paste0(colnames(test_performance_spec_youden),"-",strsplit(data,".txt")[[1]][1])

test_performance_recall_youden <- test_performance[,c("recall_youden"),drop = FALSE]
colnames(test_performance_recall_youden) <- paste0(colnames(test_performance_recall_youden),"-",strsplit(data,".txt")[[1]][1])

test_performance_precision_youden <- test_performance[,c("precision_youden"),drop = FALSE]
colnames(test_performance_precision_youden) <- paste0(colnames(test_performance_precision_youden),"-",strsplit(data,".txt")[[1]][1])

test_performance_recall_100spec <- test_performance[,c("recall_100.spec"),drop = FALSE]
colnames(test_performance_recall_100spec) <- paste0(colnames(test_performance_recall_100spec),"-",strsplit(data,".txt")[[1]][1])



all_test_probs <- cbind(all_test_probs,test_probs)
all_train_probs <- cbind(all_train_probs,train_probs)
all_test_performance_AUC <- cbind(all_test_performance_AUC,test_performance_AUC)
all_test_performance_spec_youden <- cbind(all_test_performance_spec_youden,test_performance_spec_youden)
all_test_performance_recall_youden <- cbind(all_test_performance_recall_youden,test_performance_recall_youden)
all_test_performance_precision_youden <- cbind(all_test_performance_precision_youden,test_performance_precision_youden)
all_test_performance_recall_100spec <- cbind(all_test_performance_recall_100spec,test_performance_recall_100spec)
}

write.table(all_train_probs,paste0(outdir,"/",comparison,"/",j,"_train_prob.txt"),sep = "\t", quote = FALSE ,row.names = FALSE)
write.table(all_test_probs,paste0(outdir,"/",comparison,"/",j,"_test_prob.txt"),sep = "\t", quote = FALSE ,row.names = FALSE)
write.table(all_test_performance_AUC,paste0(outdir,"/",comparison,"/",j,"_test_AUC.txt"),sep = "\t", quote = FALSE ,row.names = FALSE)
write.table(all_test_performance_spec_youden,paste0(outdir,"/",comparison,"/",j,"_test_spec_youden.txt"),sep = "\t", quote = FALSE ,row.names = FALSE)
write.table(all_test_performance_recall_youden,paste0(outdir,"/",comparison,"/",j,"_test_recall_youden.txt"),sep = "\t", quote = FALSE ,row.names = FALSE)
write.table(all_test_performance_precision_youden,paste0(outdir,"/",comparison,"/",j,"_test_precision_youden.txt"),sep = "\t", quote = FALSE ,row.names = FALSE)
write.table(all_test_performance_recall_100spec,paste0(outdir,"/",comparison,"/",j,"_test_recall_100spec.txt"),sep = "\t", quote = FALSE ,row.names = FALSE)

late_integrated <- as.data.frame(matrix(numeric(0),ncol=7))
colnames(late_integrated) <- c("comparison","seed","AUC","spec_youden","recall_youden","precision_youden","recall_100%spec")
coeff_all <- as.data.frame(matrix(numeric(0),ncol=1,nrow=length(datas)))
colnames(coeff_all) <- "Models"
coeff_all$Models <- datas
for(seed in unique(all_train_probs$seed)){
# A is data for train, B is data for test
A <- as.matrix(all_train_probs[all_train_probs$seed==seed,])
A_numeric <- matrix(as.numeric(A[,-c(1:3)]),ncol = ncol(A[,-c(1:3)]))
colnames(A_numeric) <- colnames(A[,-c(1:3)])
A_label <- as.character(A[,3])
A_label <- gsub(strsplit(comparison,"_")[[1]][1],"1",A_label)
A_label <- gsub(strsplit(comparison,"_")[[1]][3],"0",A_label)
A_label <- as.numeric(A_label)
cvfit <- cv.glmnet(A_numeric, A_label, standardize = TRUE , type.measure = "deviance", family = "binomial" , nfolds = 5, alpha = 1) 
cvfit
cvfit$lambda.min
cvfit$lambda.1se
coefficient <- coef(cvfit , s = "lambda.1se")
coeff <- as.data.frame(as.matrix(coefficient))
coeff<- coeff[-1,,drop = FALSE]
colnames(coeff) <- paste0(seed,"_",comparison)
print(coeff)
coeff_all <- cbind(coeff_all,coeff)

B <- as.matrix(all_test_probs[all_test_probs$seed==seed,])
B_numeric <- matrix(as.numeric(B[,-c(1:3)]),ncol = ncol(B[,-c(1:3)]))
colnames(B_numeric) <- colnames(B[,-c(1:3)])
B_label <- as.character(B[,3])
B_label <- gsub(strsplit(comparison,"_")[[1]][1],"1",B_label)
B_label <- gsub(strsplit(comparison,"_")[[1]][3],"0",B_label)
B_label <- as.numeric(B_label)
B_prob <- predict(cvfit,newx = B_numeric, s = "lambda.1se", type = "response")
B_prob <- as.data.frame(B_prob)
colnames(B_prob) <- 1
rownames(B_prob) <- B[,2]
#AUC
test_roc.curve <- roc(B_label,B_prob[,1],levels=c(0,1),direction=c("<"))
test_AUC <- test_roc.curve$auc[1]
print(test_AUC)
#other metric
record <- ci.coords(test_roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
spec_youden <- record$specificity[2]
recall_youden <-record$recall[2]
precision_youden <-record$precision[2]
#recall at 100% specificity
ref <- gsub("1","positive",B_label)
ref <- gsub("0","negative",ref)
cutoff = 0.5
cutoff = max(B_prob[ref=="negative",1])
i=1
pred={}
while(i<=length(B_prob[,1])){
  if(B_prob[i,1]>cutoff){
    pred[i] <- "positive"
  } else {
    pred[i] <- "negative"
  }
  i=i+1
}

pred <- factor(pred,levels = c("positive","negative"))
ref <- factor(ref,levels = c("positive","negative"))
confusion_matrix <- confusionMatrix(pred,ref)

recall_100 <- confusion_matrix$byClass[1]

late_integrated_tmp <- data.frame("comparison"=comparison,"seed"=seed,"AUC"=test_AUC,"spec_youden"=spec_youden,"recall_youden"=recall_youden,"precision_youden"=precision_youden,"recall_100%spec"=recall_100)
late_integrated <- rbind(late_integrated,late_integrated_tmp)
}

write.table(coeff_all,paste0(outdir,"/",comparison,"/",j,"_coefficients_for_different_alterations.txt"),sep = "\t", quote = FALSE ,row.names = FALSE)
write.table(late_integrated,paste0(outdir,"/",comparison,"/",j,"_late_integrated_AUC.txt"),sep = "\t", quote = FALSE ,row.names = FALSE)
}
j=j+1
}
