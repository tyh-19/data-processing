setwd("C:/Users/Tao/Desktop/")

LUAD_df <- read.csv("test.csv",sep=",",header=TRUE,encoding="UTF-8")
head(LUAD_df)

#Binomal_pvalue, the same as python function: statsmodels.stats.proportion.binom_test
#notion: need to set one group as ref, then get the same results in science,2020
#settled one group (A/B)
#group A
output_pvalue=c()
for(i in 1:nrow(LUAD_df)){
  if(as.numeric(LUAD_df[i,'Smokers'])>as.numeric(LUAD_df[i,'Never.smokers'])) {
    p <- binom.test(x=as.numeric(LUAD_df[i,'Positive.in.Never']),n=as.numeric(LUAD_df[i,'Total.Never']),p=as.numeric(LUAD_df[i,'Smokers']),alternative="less")
    output_pvalue=rbind(output_pvalue,as.data.frame(p$p.value))}
  else{
    p <- binom.test(x=as.numeric(LUAD_df[i,'Positive.in.Never']),n=as.numeric(LUAD_df[i,'Total.Never']),p=as.numeric(LUAD_df[i,'Smokers']),alternative="greater")
    output_pvalue=rbind(output_pvalue,as.data.frame(p$p.value))
  }
}
rownames(output_pvalue) <- LUAD_df$﻿taxo
colnames(output_pvalue) <- "2nd-Binomal-p-A"
output_pvalue_A <- cbind(LUAD_df,output_pvalue)
#group B
output_pvalue=c()
for(i in 1:nrow(LUAD_df)){
  if(as.numeric(LUAD_df[i,'Smokers'])>as.numeric(LUAD_df[i,'Never.smokers'])) {
    p <- binom.test(x=as.numeric(LUAD_df[i,'Positive.in.Smokers']),n=as.numeric(LUAD_df[i,'Total.Smokers']),p=as.numeric(LUAD_df[i,'Never.smokers']),alternative="greater")
    output_pvalue=rbind(output_pvalue,as.data.frame(p$p.value))}
  else{
    p <- binom.test(x=as.numeric(LUAD_df[i,'Positive.in.Smokers']),n=as.numeric(LUAD_df[i,'Total.Smokers']),p=as.numeric(LUAD_df[i,'Never.smokers']),alternative="less")
    output_pvalue=rbind(output_pvalue,as.data.frame(p$p.value))
  }
}
rownames(output_pvalue) <- LUAD_df$﻿taxo
colnames(output_pvalue) <- "2nd-Binomal-p-B"
output_pvalue <- cbind(output_pvalue_A,output_pvalue)
#write.csv(output_pvalue,"test.csv")





###############################freezed##########################################################
################################################################################################

#2 proportion_pvalue
#different from python function: statsmodels.stats.proportion.proportions_ztest

output_pvalue_2prop=c()
for(i in 1:nrow(LUAD_df)){
  p <- prop.test(c(as.numeric(LUAD_df[i,'Positive.in.Smokers']),as.numeric(LUAD_df[i,'Positive.in.Never'])),c(as.numeric(LUAD_df[i,'Total.Smokers']),as.numeric(LUAD_df[i,'Total.Never'])))
  output_pvalue_2prop=rbind(output_pvalue_2prop,as.data.frame(p$p.value))
}
rownames(output_pvalue_2prop) <- LUAD_df$﻿taxo
colnames(output_pvalue_2prop) <- "2nd-2-proportion-p"
head(output_pvalue_2prop)
output_pvalue_2prop <- cbind(output_pvalue,output_pvalue_2prop)
#write.csv(output_pvalue_2prop,"p_2prop.csv")

#phyper test(hypergeometry test)
#phyper(q,M,N-M,k)
#q=picked number of a specific genus(N=genus exist in one group)
#M=exist number of the genus in both group(M=genus exist in group A+B)
#N=all samples
#k=trails number(k>q,k=sample number in one group)
N=as.numeric(LUAD_df[i,'Total.Smokers'])+as.numeric(LUAD_df[i,'Total.Never'])
output_pvalue_phyper=c()
for(i in 1:nrow(LUAD_df)) {
  if(as.numeric(LUAD_df[i,'Smokers'])>as.numeric(LUAD_df[i,'Never.smokers'])) {
    M=as.numeric(LUAD_df[i,'Positive.in.Smokers'])+as.numeric(LUAD_df[i,'Positive.in.Never'])
    p <- phyper(as.numeric(LUAD_df[i,'Positive.in.Smokers'])-1,M,N-M,as.numeric(LUAD_df[i,'Total.Smokers']),lower.tail=FALSE)
    output_pvalue_phyper=rbind(output_pvalue_phyper,as.data.frame(p))
    }
  else{
    M=as.numeric(LUAD_df[i,'Positive.in.Smokers'])+as.numeric(LUAD_df[i,'Positive.in.Never'])
    p <- phyper(as.numeric(LUAD_df[i,'Positive.in.Never'])-1,M,N-M,as.numeric(LUAD_df[i,'Total.Never']),lower.tail=FALSE)
    output_pvalue_phyper=rbind(output_pvalue_phyper,as.data.frame(p))
    }
}
rownames(output_pvalue_phyper) <- LUAD_df$﻿taxo
colnames(output_pvalue_phyper) <- "2nd-hypergeometry-p"
head(output_pvalue_phyper)
output_pvalue_phyper <- cbind(output_pvalue_2prop,output_pvalue_phyper)
write.csv(output_pvalue_phyper,"test.csv")

#fdr
pvalue_to_fdr <- read.csv("LUAD_pvalue_filtered.csv",header=T)
pvalue_to_fdr <- as.numeric(unlist(pvalue_to_fdr$p.p.value))
length(pvalue_to_fdr)

fdr <- p.adjust(pvalue_to_fdr, method = "BH", n = length(pvalue_to_fdr))
as.data.frame(fdr)
pvalue_to_fdr <- cbind(pvalue_to_fdr,fdr)
write.csv(pvalue_to_fdr,"LUAD_filtered.csv")

??read.csv
length(pvalue_to_fdr$p.p.value)

##test
i=1
M=as.numeric(LUAD_df[i,'Postive.in.Current.smoker'])+as.numeric(LUAD_df[i,'Postive.in.Never.smoker'])
p <- phyper(as.numeric(LUAD_df[i,'Postive.in.Current.smoker'])-1,M,34-M,as.numeric(LUAD_df[i,'Total.Current.smoker']),lower.tail=FALSE)
