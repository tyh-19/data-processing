#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))


parser <- ArgumentParser(description='Agregate rMATs output stats and inclevel for A5SS,A3SS,MXE,RI,SE')
parser$add_argument('-i', '--input_dir', type='character', required=TRUE,
    help='input directory, which contains stats and inclevel for A5SS,A3SS,MXE,RI,SE(these files generate from yunfans script)')
parser$add_argument('-m', '--output_matrix', type='character', required=TRUE,
    help='output inclevel matrix.')
parser$add_argument('-s', '--output_stat', type='character', required=TRUE,
    help='output differential splicing stats.')
args <- parser$parse_args()

    file_tosum <- args$input_dir
    files <- dir(file_tosum)
    stats <- grep("stats",files,value = TRUE)
    inclevels <- grep("inc_level",files,value = TRUE)
    output_stat <- args$output_stat
    output_inclevel <- args$output_matrix
    
    stat_all <- data.frame()
    for(i in stats){
      stat <- read.csv(paste0(file_tosum,"/",i),sep = "\t", header = TRUE, row.names = 1)
      stat_all <- rbind(stat_all,stat)
    }
    write.table(stat_all,output_stat,quote = FALSE, sep = "\t")
    
    inclevel_matrix <- data.frame()
    for(i in inclevels){
      inclevel <- read.csv(paste0(file_tosum,"/",i),sep = "\t", header = TRUE, row.names = 1)
      inclevel_matrix <- rbind(inclevel_matrix,inclevel)
    }
    write.table(inclevel_matrix,output_inclevel,quote = FALSE, sep = "\t")
