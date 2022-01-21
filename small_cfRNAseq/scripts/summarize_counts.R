#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))


parser <- ArgumentParser(description='Summarize count')
parser$add_argument('-i', '--input_dir', type='character', required=TRUE,
    help='input directory, which contains count files(end with .txt)')
parser$add_argument('-t', '--type', type='character', required=TRUE,
    help='data type.')
parser$add_argument('-o', '--output_matrix', type='character', required=TRUE,
    help='output matrix.')
args <- parser$parse_args()

  
  input_dir <- args$input_dir
  output_matrix <- args$output_matrix
  files <- dir(input_dir)
  type <- args$type
  
  i=1
  count_matrix <- as.data.frame(matrix(numeric(0),ncol=1))
  colnames(count_matrix) <- c("feature") 
  count_matrix$`feature` <- as.factor(count_matrix$`feature`)
  while(i<=length(files)){
  count <- read.table(paste0(input_dir,"/",files[i],"/",type,".txt"))
  sample <- as.character(lapply(strsplit(files[i],"_1",fixed = TRUE),function(x) x[1]))
  colnames(count) <- c("feature",sample)
  count_matrix <- full_join(count_matrix,count,by=c("feature"="feature"))
  i=i+1
  }
  
  rownames(count_matrix) <- count_matrix$feature
  count_matrix <- count_matrix[,-which(colnames(count_matrix)=="feature")]
  write.table(count_matrix,paste0(output_matrix),sep = "\t",quote = FALSE)
