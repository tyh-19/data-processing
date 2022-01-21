#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='Summarize matrics')
parser$add_argument('-i', '--input_dir', type='character', required=TRUE,
    help='input dir')
parser$add_argument('-p', '--positive', type='character', required=TRUE,
    choices=c('CRC','STAD','NC'))
parser$add_argument('-n','--negative', type='character', required=TRUE,
    choices=c('CRC','STAD','NC'))
parser$add_argument('-o', '--output_dir', type='character', required=TRUE,
    help='output dir')
args <- parser$parse_args()

input_dir <- args$input_dir
  output_dir_initial <- args$output_dir
  positive <- args$positive
  negative <- args$negative
  
  
  full_alterations <- dir(input_dir)
  full_alterations <- full_alterations[-grep("Merged",full_alterations)]
  #full_alterations <- full_alterations[-grep("txt",full_alterations)]
  full_alterations <- full_alterations[full_alterations!="CNV"]
  full_alterations <- full_alterations[full_alterations!="Splicing_FDR0.05"]
  
  n=1
  report <- data.frame(matrix(numeric(0),ncol = 7, nrow = length(full_alterations)))
  while(n <= length(full_alterations)){
  colnames(report) <- c("Alteration","Input_dir","AUC","Specificity_0.95specificity","Recall_0.95specificity","Specificity_youden","Recall_youden")
  AUC <- read.csv(paste0(input_dir,"/",full_alterations[n],"/",positive,"vs",negative,"/AUC.txt"),sep = "\t")
  recall_0.95 <- read.csv(paste0(input_dir,"/",full_alterations[n],"/",positive,"vs",negative,"/Sensitivity_0.95specificity.txt"),sep = "\t")
  specificity_0.95 <- read.csv(paste0(input_dir,"/",full_alterations[n],"/",positive,"vs",negative,"/Specificity_0.95specificity.txt"),sep = "\t")
  recall_youden <- read.csv(paste0(input_dir,"/",full_alterations[n],"/",positive,"vs",negative,"/Sensitivity_youden.txt"),sep = "\t")
  specificity_youden <- read.csv(paste0(input_dir,"/",full_alterations[n],"/",positive,"vs",negative,"/Sensitivity_youden.txt"),sep = "\t")
  report[n,which(colnames(report)=="Alteration")] <- full_alterations[n]
  report[n,which(colnames(report)=="Input_dir")] <- input_dir
  report[n,which(colnames(report)=="AUC")] <- AUC[1,1]
  report[n,which(colnames(report)=="Specificity_0.95specificity")] <- specificity_0.95[1,1]
  report[n,which(colnames(report)=="Recall_0.95specificity")] <- recall_0.95[1,1]
  report[n,which(colnames(report)=="Specificity_youden")] <- specificity_youden[1,1]
  report[n,which(colnames(report)=="Recall_youden")] <- recall_youden[1,1]
  n=n+1
  }
  
  write.table(report,paste0(output_dir_initial,"/report.txt"),quote = FALSE,sep = "\t")

