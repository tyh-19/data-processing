install.packages("kknn")
library(kknn)
# train test data_for_prediction
QC_raw <- read.csv("C:/Users/Tao/Desktop/远程办公文件/lulab 2 (1)/lulab/Liquid Biopsy/02.qPCR QC/QC_raw_data.csv",header=TRUE)
head(QC_raw)
QC_raw$quality
QC_test
QC_test_trueresult
QC_train
QC_pengfei_for_prediction

#model
KNN_model_4 <- kknn(quality~miR15+miR16+miR21+miR451A+miR468+hY4+miR1228+actb+gapdh+HULC.3+MALAT1,QC_train, QC_test, k= 4, distance = 2, kernel = "optimal", ykernel = NULL, scale=TRUE, contrasts= c('unordered' = "contr.dummy", ordered ="contr.ordinal"))
KNN_model_4 <- kknn(quality~miR15+miR16+miR21+miR451A+miR468+hY4+miR1228,QC_raw, QC_pengfei_for_prediction, k= 4, distance = 2, kernel = "optimal", ykernel = NULL, scale=TRUE, contrasts= c('unordered' = "contr.dummy", ordered ="contr.ordinal"))
KNN_model_4 <- kknn(quality~miR15+miR16+miR21+actb+gapdh,QC_train, QC_test, k= 4, distance = 2, kernel = "optimal", ykernel = NULL, scale=TRUE, contrasts= c('unordered' = "contr.dummy", ordered ="contr.ordinal"))
#evaluation
KNN_model_4
fit <- fitted(KNN_model_4) 
fit
QC_pengfei_for_prediction$quality_KNN <- fit
QC_pengfei_for_prediction
table(fit,QC_pengfei_for_prediction$quality)
write.csv(QC_pengfei_for_prediction,"C:/Users/Tao/Desktop/远程办公文件/lulab 2 (1)/lulab/Liquid Biopsy/02.qPCR QC/prediction/QC_pengfei_predicted.csv")

class(KNN_model_4$prob)
prediction_KNN <- data.frame(KNN_model_4$prob)
prediction_KNN
#ROC
library(ROCR)
pred <- prediction(prediction_KNN$good,QC_test_trueresult$quality)
perf <- performance(pred,"tpr","fpr")
plot(perf, col='blue',lty=2)
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

#PCA feature_selection
##princomp
QC_raw[,2:12]
QC_raw_PCA <- princomp(QC_raw[,2:12],cor=TRUE)
summary(QC_raw_PCA,loadings=TRUE)
predict(QC_raw_PCA)
screeplot(QC_raw_PCA,type="lines")
biplot(QC_raw_PCA,choices=1:2,scale=1,pc.biplot=FALSE)
scoresdata=QC_raw_PCA$scores
head(scoresdata)
##prcomp
QC_raw_PCA_ggplot <- prcomp(QC_raw[,2:12], scale = TRUE)
head(QC_raw_PCA_ggplot,1)
head(QC_raw_PCA,1)
##princomp and prcomp are the same
class(QC_raw_PCA_ggplot)
QC_raw_PCA$scores
QC_raw_pcs <- data.frame(QC_raw_PCA$scores,quality=QC_raw$quality)

head(QC_raw_pcs)
write.csv(QC_raw_pcs,"C:/Users/Tao/Desktop/远程办公文件/lulab 2 (1)/lulab/Liquid Biopsy/02.qPCR QC/prediction/QC_raw_pcs.csv")

##percentage
percentage<-round(QC_raw_PCA$sdev/sum(QC_raw_PCA$sdev) * 100,2)
percentage<-paste(colnames(QC_raw_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

##ggplot
ggplot(QC_raw_pcs,aes(x=Comp.1,y=Comp.2,z=Comp.3))+
  geom_point(aes(colour=quality))+
  xlab(percentage[1])+ylab(percentage[2])+
  theme_bw()
  #theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line= element_line(colour = "black"))+
  #facet_wrap(~Operator)

#annotate('text', label = '14', x = -2, y = -1.25, size = 5, colour = '#f8766d')
#stat_ellipse(level = 0.95)+

##contribution of each feature
QC_raw_PCA_ggplot$rotation
PCA_r <- as.data.frame(QC_raw_PCA_ggplot$rotation)
PCA_r
PCA_r$feature <- row.names(PCA_r)
ggplot(PCA_r,aes(x=PC1,y=PC2,label=feature,color=feature )) + geom_point()+ geom_text(size=3,hjust=0.6,vjust=-0.7) + theme_classic()

