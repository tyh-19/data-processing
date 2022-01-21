 #! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(doParallel))
parser <- ArgumentParser(description='Summarize SNP allele frequency from vcf.')
parser$add_argument('-a', '--input_Alt_matrix', type='character', required=TRUE,
    help='input Alt count')
parser$add_argument('-r', '--input_Ref_matrix', type='character', required=TRUE,
    help='input Ref count')
parser$add_argument('-p', '--positive', type='character', required=TRUE,
    help='input positive group prefix')
parser$add_argument('-n', '--negative', type='character', required=TRUE,
    help='input negative group prefix')
parser$add_argument('-o', '--output_table', type='character', required=TRUE,
    help='output fisher exact test result')
args <- parser$parse_args()

  cl <- makeCluster(16) #not to overload your computer
  registerDoParallel(cl)
  
  final_SNP_Alt <- read.csv(args$input_Alt_matrix,sep = "\t",header = TRUE,row.names = 1,stringsAsFactors = FALSE)
  colnames(final_SNP_Alt) <- gsub(".","-",fixed=TRUE,colnames(final_SNP_Alt))
  colnames(final_SNP_Alt) <- gsub("_Alt","",colnames(final_SNP_Alt))
  final_SNP_Ref <- read.csv(args$input_Ref_matrix,sep = "\t",header = TRUE,row.names = 1,stringsAsFactors = FALSE)
  colnames(final_SNP_Ref) <- gsub(".","-",fixed=TRUE,colnames(final_SNP_Ref))
  colnames(final_SNP_Ref) <- gsub("_Ref","",colnames(final_SNP_Ref))

  positive <- args$positive
  negative <- args$negative
  
  positive_id <- grep(positive,colnames(final_SNP_Alt))
  negative_id <- grep(negative,colnames(final_SNP_Ref))
  
  message('Positive sample number:', length(positive_id))
  message('Negative sample number:', length(negative_id))

  SNP_Alt_count <- data.frame("SNP"=rownames(final_SNP_Alt),
                              "Positive_Alt_count"=rowSums(final_SNP_Alt[,positive_id]),
                              "Negative_Alt_count"=rowSums(final_SNP_Alt[,negative_id]))
  
  SNP_Ref_count <- data.frame("SNP"=rownames(final_SNP_Ref),
                              "Positive_Ref_count"=rowSums(final_SNP_Ref[,positive_id]),
                              "Negative_Ref_count"=rowSums(final_SNP_Ref[,negative_id]))
  j=1
  SNP_DE <- data.frame(matrix(numeric(0),ncol=5,nrow=nrow(SNP_Ref_count)))
  colnames(SNP_DE) <- c("SNP","Positive group AF","Negative group AF","pvalue","padj")
  
  total <- nrow(SNP_Ref_count)
  while(j<=nrow(SNP_Ref_count)) {
  #get SNP_site
  SNP_DE[j,"SNP"] <- as.character(SNP_Alt_count[j,1])
  
  #assemble two-dimensional contingency table
  Alt_count <- SNP_Alt_count[j,-1]
  colnames(Alt_count) <- c(positive,negative)
  rownames(Alt_count) <- "Alt"
  Ref_count <- SNP_Ref_count[j,-1]
  colnames(Ref_count) <- c(positive,negative)
  rownames(Ref_count) <- "Ref"
  
  SNP_fishertest <- as.matrix(rbind(Alt_count,Ref_count))
  fisher_result <- fisher.test(SNP_fishertest)
  SNP_DE[j,"pvalue"] <- fisher_result$p.value
  
  AF <- Alt_count/(Ref_count+Alt_count)
  SNP_DE[j,"Positive group AF"] <- AF[1,positive]
  SNP_DE[j,"Negative group AF"] <- AF[1,negative]
  message('Progress:',j,'/',total)
  j=j+1
  }
  
  SNP_DE$padj <- p.adjust(SNP_DE$pvalue,method = "BH")
  
  rownames(SNP_DE) <- SNP_DE$SNP
  
  write.table(SNP_DE[,-which(colnames(SNP_DE)=="SNP")],args$output_table,sep = "\t",quote = FALSE)

  stopCluster(cl)  
