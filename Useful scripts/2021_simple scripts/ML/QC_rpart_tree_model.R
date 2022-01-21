QC_raw <- read.csv("C:/Users/Tao/Desktop/远程办公文件/lulab 2 (1)/lulab/Liquid Biopsy/02.qPCR QC/QC_raw_data.csv",header=TRUE)
head(QC_raw)
install.packages("rpart")
library(rpart)
row.names(QC_raw) <- QC_raw$锘縄D
row.names(QC_raw)
QC_test <- QC_raw[1:5,1:12]
QC_test_trueresult <- QC_raw[1:5,] 
QC_train <- QC_raw[6:26,]
QC_train
tree_model <- rpart(quality~miR15+miR16+miR21+miR451A+miR468+hY4+miR1228+actb+gapdh+HULC.3+MALAT1,method="class",
             data=QC_train)
plot(tree_model,uniform=TRUE,main="Classification Tree for QC_raw[6:26,]")
tree_model_text <- text(tree_model,use.n=TRUE,all=TRUE)
tree_model_text


qPCR_quality_prediected <- predict(tree_model,QC_test,type="class") 
QC_test
qPCR_quality_prediected
QC_test_trueresult


##test QC_pengfei_for_prediction
QC_pengfei_for_prediction <- read.csv("C:/Users/Tao/Desktop/远程办公文件/lulab 2 (1)/lulab/Liquid Biopsy/02.qPCR QC/prediction/QC_pengfei_for prediction.csv",header=TRUE)
head(QC_pengfei_for_prediction)
QC_pengfei_for_prediction_predicted <- predict(tree_model,QC_pengfei_for_prediction,type="class") 
QC_pengfei_for_prediction_predicted
QC_pengfei_for_prediction$quality <- QC_pengfei_for_prediction_predicted
head(QC_pengfei_for_prediction)
colnames(QC_pengfei_for_prediction)
class(QC_pengfei_for_prediction$quality)
QC_pengfei_for_prediction$quality <- as.factor(QC_pengfei_for_prediction$quality)
write.csv(QC_pengfei_for_prediction,"C:/Users/Tao/Desktop/远程办公文件/lulab 2 (1)/lulab/Liquid Biopsy/02.qPCR QC/prediction/QC_pengfei_predicted.csv")
