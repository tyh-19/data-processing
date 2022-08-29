#library
library(edgeR)
library(dplyr)
library(progress)
#library()
#function
Wilcox_test <- function(mat_raw,positive_samples,negative_samples){
  norm_method <- 'NA'
  positive <- as.character(positive_samples)
  negative <- as.character(negative_samples)
  group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
  
  if(norm_method == 'NA' ){
    message('Matrix output without normalization.')
    mat_raw[is.na(mat_raw)] <- 0
    matrix <- mat_raw
  }else{
    message('Matrix output normalized by:',norm_method)
    matrix <- cpm(mat_raw, method=norm_method)
  }
  
  colnames(matrix) <- gsub(".","-",fixed = TRUE,colnames(matrix))
  matrix <- matrix[,as.character(c(positive_samples,negative_samples))]
  
  
  
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
  deltaAF <- apply(matrix[,positive], 1, mean) -
    apply(matrix[,negative], 1, mean)
  
  positive_gini <- as.numeric(gini(t(matrix[,positive])))
  negative_gini <- as.numeric(gini(t(matrix[,negative])))
  
  res <- data.frame(log2FoldChange=logFC,
                    deltaAF=deltaAF,
                    pvalue=pvalues, 
                    padj=p.adjust(pvalues, method='BH'),
                    baseMean=apply(matrix, 1, mean),
                    positive_gini = positive_gini,
                    negative_gini = negative_gini)
  #write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE)
  res
}

#binary feature
{
#select complementary features
alteration1 <- "cfRNA abundance"
alteration2 <- "cfDNA methylation"

output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220713_yuhuan/complementary_events"
dir.create(output_dir,recursive = TRUE)
data1 <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/Plasma_TPM.txt",sep = "\t",check.names = FALSE,header = TRUE, row.names = 1)
data2 <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/DNA/CPM-TMM_matrix_promoter.txt",sep = "\t",check.names = FALSE,header = TRUE, row.names = 1)
colnames(data1) <- gsub("-pico","",colnames(data1))
colnames(data2) <- gsub("-me","",colnames(data2))

samples <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/ML_samples.csv",row.names = NULL)
positive <- "CRC|STAD"
negative <- "HD"
positive_samples <- as.character(samples[grep(positive,samples$Group),]$sample_id)
n_positive <- length(positive_samples)
negative_samples <- as.character(samples[grep(negative,samples$Group),]$sample_id)
n_negative <- length(negative_samples)

data1_diff <- Wilcox_test(data1,positive_samples,negative_samples)
data1_diff <- na.omit(data1_diff)
data1_selected <- rownames(data1_diff[data1_diff$pvalue<0.000001,])
data2_diff <- Wilcox_test(data2,positive_samples,negative_samples)
data2_diff <- na.omit(data2_diff)
data2_selected <- rownames(data2_diff[data2_diff$pvalue<0.001,])


comb_events <- expand.grid(data1_selected,data2_selected)
message(nrow(comb_events))
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = nrow(comb_events), clear = FALSE, width= 60)

i=1
result_final <- as.data.frame(matrix(numeric(0),ncol=16))
colnames(result_final) <- c("data1","data2","event1","event2","cutoff_1_upper","cutoff_1_lower","cutoff_2_upper","cutoff_2_lower",
                            "Upper_upper_union","Upper_upper_intersect","Upper_lower_union","Upper_lower_intersect",
                            "Lower_upper_union","Lower_upper_intersect","Lower_lower_union","Lower_lower_intersect")
result <- as.data.frame(matrix(numeric(0),ncol=16))
colnames(result) <- c("data1","data2","event1","event2","cutoff_1_upper","cutoff_1_lower","cutoff_2_upper","cutoff_2_lower",
                            "Upper_upper_union","Upper_upper_intersect","Upper_lower_union","Upper_lower_intersect",
                            "Lower_upper_union","Lower_upper_intersect","Lower_lower_union","Lower_lower_intersect")
while(i<=nrow(comb_events)){
#while(i<=1000){
event1 <- as.character(comb_events[i,1])
positive_data1 <- as.data.frame(t(data1[event1,as.character(positive_samples)]))
negative_data1 <- as.data.frame(t(data1[event1,as.character(negative_samples)]))
cutoff_1_upper <- negative_data1[order(negative_data1[,event1],decreasing = TRUE),][2]
cutoff_1_lower <- negative_data1[order(negative_data1[,event1],decreasing = TRUE),][n_negative-1]

event2 <- as.character(comb_events[i,2])
positive_data2 <- as.data.frame(t(data2[event2,as.character(positive_samples)]))
negative_data2 <- as.data.frame(t(data2[event2,as.character(negative_samples)]))
cutoff_2_upper <- negative_data2[order(negative_data2[,event2],decreasing = TRUE),][2]
cutoff_2_lower <- negative_data2[order(negative_data2[,event2],decreasing = TRUE),][n_negative-1]

data1_outlier_upper <- rownames(positive_data1[positive_data1[[event1]]>cutoff_1_upper,,drop=FALSE])
data1_outlier_lower <- rownames(positive_data1[positive_data1[[event1]]<cutoff_1_lower,,drop=FALSE])
data2_outlier_upper <- rownames(positive_data2[positive_data2[[event2]]>cutoff_2_upper,,drop=FALSE])
data2_outlier_lower <- rownames(positive_data2[positive_data2[[event2]]<cutoff_2_lower,,drop=FALSE])

Upper_upper_union <- length(union(data1_outlier_upper,data2_outlier_upper))/n_positive
Upper_upper_intersect <- length(intersect(data1_outlier_upper,data2_outlier_upper))/n_positive

Upper_lower_union <- length(union(data1_outlier_upper,data2_outlier_lower))/n_positive
Upper_lower_intersect <- length(intersect(data1_outlier_upper,data2_outlier_lower))/n_positive

Lower_upper_union <- length(union(data1_outlier_lower,data2_outlier_upper))/n_positive
Lower_upper_intersect <- length(intersect(data1_outlier_lower,data2_outlier_upper))/n_positive
  
Lower_lower_union <- length(union(data1_outlier_lower,data2_outlier_lower))/n_positive
Lower_lower_intersect <- length(intersect(data1_outlier_lower,data2_outlier_lower))/n_positive

result_tmp <- data.frame("data1"=alteration1,"data2"=alteration2,"event1"=event1,"event2"=event2,"cutoff_1_upper"=cutoff_1_upper,"cutoff_1_lower"=cutoff_1_lower,"cutoff_2_upper"=cutoff_2_upper,"cutoff_2_lower"=cutoff_2_lower,
                         "Upper_upper_union"=Upper_upper_union,"Upper_upper_intersect"=Upper_upper_intersect,"Upper_lower_union"=Upper_lower_union,"Upper_lower_intersect"=Upper_lower_intersect,
                         "Lower_upper_union"=Lower_upper_union,"Lower_upper_intersect"=Lower_upper_intersect,"Lower_lower_union"=Lower_lower_union,"Lower_lower_intersect"=Lower_lower_intersect)
result <- rbind(result,result_tmp)

if(i%%5000==0){
  result_final <- rbind(result_final,result)
  write.csv(result,paste0(output_dir,"/",i/5000,"_Result5000.csv"))
  result <- as.data.frame(matrix(numeric(0),ncol=16))
  colnames(result) <- c("data1","data2","event1","event2","cutoff_1_upper","cutoff_1_lower","cutoff_2_upper","cutoff_2_lower",
                              "Upper_upper_union","Upper_upper_intersect","Upper_lower_union","Upper_lower_intersect",
                              "Lower_upper_union","Lower_upper_intersect","Lower_lower_union","Lower_lower_intersect")
  i=i+1
} else if(i==(nrow(comb_events)-1)){
  result_final <- rbind(result_final,result)
  write.csv(result,paste0(output_dir,"/",i%/%5000+1,"_Result",i%%5000,".csv"))
  i=i+1
} else {
  i=i+1
}
pb$tick()
Sys.sleep(1 / 100)
}
  
write.csv(result_final,paste0(output_dir,"/",positive,"_vs_",negative,"|",alteration1,"|",alteration2,"|complementary_events.csv"))
}

#single feature
{
  #select complementary features
  alteration3 <- "cfDNA copy number"
  
  output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220713_yuhuan/complementary_events"
  dir.create(output_dir,recursive = TRUE)
  data3 <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/DNA/CPM-TMM_matrix_gene.correctGC.txt",sep = "\t",check.names = FALSE,header = TRUE, row.names = 1)
  colnames(data3) <- gsub("-wgs","",colnames(data3))
  
  samples <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/ML_samples.csv",row.names = NULL)
  positive <- "STAD"
  negative <- "CRC"
  positive_samples <- as.character(samples[grep(positive,samples$Group),]$sample_id)
  n_positive <- length(positive_samples)
  negative_samples <- as.character(samples[grep(negative,samples$Group),]$sample_id)
  n_negative <- length(negative_samples)
  
  data3_diff <- Wilcox_test(data3,positive_samples,negative_samples)
  data3_diff <- na.omit(data3_diff)
  data3_selected <- rownames(data3_diff[data3_diff$pvalue<0.00001,])
  
  
  message(length(data1_selected))
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = length(data3_selected)-1, clear = FALSE, width= 60)
  
  i=1
  result_final <- as.data.frame(matrix(numeric(0),ncol=6))
  colnames(result_final) <- c("data3","event3","cutoff_3_upper","cutoff_3_lower",
                              "Upper","Lower")
  result <- as.data.frame(matrix(numeric(0),ncol=6))
  colnames(result) <- c("data3","event3","cutoff_3_upper","cutoff_3_lower",
                        "Upper","Lower")
  while(i<=length(data3_selected)-1){
    #while(i<=1000){
    event3 <- data3_selected[i]
    positive_data3 <- as.data.frame(t(data3[event3,as.character(positive_samples)]))
    negative_data3 <- as.data.frame(t(data3[event3,as.character(negative_samples)]))
    cutoff_3_upper <- negative_data3[order(negative_data3[,event3],decreasing = TRUE),][2]
    cutoff_3_lower <- negative_data3[order(negative_data3[,event3],decreasing = TRUE),][n_negative-1]
    
    data3_outlier_upper <- rownames(positive_data3[positive_data3[[event3]]>cutoff_3_upper,,drop=FALSE])
    data3_outlier_lower <- rownames(positive_data3[positive_data3[[event3]]<cutoff_3_lower,,drop=FALSE])
    
    Upper <- length(data3_outlier_upper)/n_positive
    
    Lower <- length(data3_outlier_lower)/n_positive
    
    result_tmp <- data.frame("data3"=alteration3,"event3"=event3,"cutoff_3_upper"=cutoff_3_upper,"cutoff_3_lower"=cutoff_3_lower,
                             "Upper"=Upper,"Lower"=Lower)
    result <- rbind(result,result_tmp)
    
    if(i%%5000==0){
      result_final <- rbind(result_final,result)
      write.csv(result,paste0(output_dir,"/",i/5000,"_Result5000.csv"))
      result <- as.data.frame(matrix(numeric(0),ncol=16))
      colnames(result) <- c("data3","event3","cutoff_3_upper","cutoff_3_lower",
                            "Upper","Lower")
      i=i+1
    } else if(i==(length(data3_selected)-1)){
      result_final <- rbind(result_final,result)
      write.csv(result,paste0(output_dir,"/",i%/%5000+1,"_Result",i%%5000,".csv"))
      i=i+1
    } else {
      i=i+1
    }
    pb$tick()
    Sys.sleep(1 / 100)
  }
  
  write.csv(result_final,paste0(output_dir,"/",positive,"_vs_",negative,"|",alteration3,"|single_events.csv"))
}

library(scatterplot3d)

#Group_color <- c(CRC=alpha("#FCB514",alpha = 1),STAD=alpha("red",alpha = 1),HD=alpha("blue",alpha = 1),GIC="#FCB514")
#event_plot <- t(rbind(data1["ENSG00000112419.14|11224",c(positive_samples,negative_samples)],data2["ENSG00000072818.11|7421",c(positive_samples,negative_samples)]))
#event_plot <- as.data.frame(event_plot)
#event_plot$group <- as.character(lapply(strsplit(rownames(event_plot),"-",fixed = TRUE), function(x) x[1]))
#event_plot$group <- gsub("NC","HD",event_plot$group)
#ggplot(event_plot,aes(x=event_plot$`ENSG00000112419.14|11224`,y=event_plot$`ENSG00000072818.11|7421`))+
#  scale_color_manual(values = Group_color)+
#  geom_vline(xintercept = 32.8624, linetype='dashed')+
#  geom_hline(yintercept = 5.2275, linetype='dashed')+
#  geom_point(aes(color = event_plot$group))+
#  theme_bw()

result_binary <- read.csv((paste0(output_dir,"/","CRC|STAD","_vs_","HD","|",alteration1,"|",alteration2,"|complementary_events.csv")))
result_single <- read.csv(paste0(output_dir,"/","STAD","_vs_","CRC","|",alteration3,"|single_events.csv"))

result_binary$Upper_upper <- result_binary$Upper_upper_union-result_binary$Upper_upper_intersect

result_binary$Upper_lower <- result_binary$Upper_lower_union-result_binary$Upper_lower_intersect

result_binary$Lower_lower <- result_binary$Lower_lower_union-result_binary$Lower_lower_intersect

result_binary$Lower_upper <- result_binary$Lower_upper_union-result_binary$Lower_upper_intersect

selected <- rbind(result_binary[order(result_binary$Upper_upper,decreasing = TRUE)[1],],
result_binary[order(result_binary$Upper_lower,decreasing = TRUE)[1],],
result_binary[order(result_binary$Lower_lower,decreasing = TRUE)[1],],
result_binary[order(result_binary$Lower_upper,decreasing = TRUE)[1],])
selected$max <- NA
selected[1,"max"] <- max(selected$Upper_upper)
selected[2,"max"] <- max(selected$Upper_lower)
selected[3,"max"] <- max(selected$Lower_lower)
selected[4,"max"] <- max(selected$Lower_upper)

event1 <- as.character(selected[order(selected$max,decreasing = TRUE)[1],]$event1)
event2 <- as.character(selected[order(selected$max,decreasing = TRUE)[1],]$event2)

event3 <- as.character(result_single[order(result_single$Upper,decreasing = TRUE)[1],]$event3)


combined_3d <- as.data.frame(t(rbind(
        data1[event1,as.character(samples$sample_id)],
        data2[event2,as.character(samples$sample_id)],
        data3[event3,as.character(samples$sample_id)]
        )))
combined_3d$group <- as.character(samples$Group)
colors <- c("#FCB514", "blue", "red")
colors <- colors[as.numeric(as.factor(combined_3d$group))]
s3d <- scatterplot3d(combined_3d[,1:3], pch=16, color = colors,angle = 90,#type="h",
                     grid=TRUE, box=TRUE)
legend("topright", legend = levels(as.factor(combined_3d$group)),
       col =  c("#FCB514", "blue", "red"), pch = 16, bg = "transparent")
text(s3d$xyz.convert(combined_3d[,1:3]), labels = rownames(combined_3d),
     cex= 0.5, col = "black")


