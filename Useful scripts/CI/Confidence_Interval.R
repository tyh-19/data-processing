library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type="character", 
              help="Input file names 'miso_Insert_length_summary.txt' or other file contain a string of numbers each column, split by tab"),
  make_option(c("-o", "--output"), type="character", 
              help="Output directory and file name"),
  make_option(c("-b", "--bootstrap_num"), type="integer", default=1000,
              help="Bootstrap number, default is 1000"),
  make_option(c("-c", "--confidence_level"), type="double", default=0.95,
              help="Confidence level, defaut is 0.95")
)

opt <- parse_args(OptionParser(option_list=option_list))

#Confidence Interval function
#x is a list of numbers
#B is sample times
#I is confidence interval, usually 0.95
#name is colnames of input length summary
Bootstrap_CI <- function(x,B,I,name){
  bstrap <- c()
  for (i in 1:B){
    bstrap <- c(bstrap,mean(sample(x,length(x),replace=T)))
  }
  mean <- mean(x)
  median <- median(x)
  sd <- sd(x)
  output <- cbind(as.data.frame(quantile(bstrap,(1-I)/2)),as.data.frame(quantile(bstrap,1-(1-I)/2)),mean-sd,mean+sd,mean,median)
  colnames(output) <-  c(paste0((1-I)/2*100,"%"),paste0((1-(1-I)/2)*100,"%"),"mean-sd","mean+sd","mean","median")
  rownames(output) <- name
  output
}

length <- read.csv(opt$input,sep = "\t",header = TRUE)
#length <- length[-1,]

i=1
all_sample_length={}

while(i<=length(colnames(length))){
sample <- na.omit(length[,i])

output <- Bootstrap_CI(sample,opt$bootstrap_num,opt$confidence_level,colnames(length)[i])

all_sample_length <- rbind(all_sample_length,output)

i=i+1
}

write.table(all_sample_length,opt$output, sep = "\t")
