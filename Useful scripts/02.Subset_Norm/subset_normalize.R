#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='Differential expression')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('-s', '--sample_ids', type='character', required=TRUE,
    help='file only contains ids, each line is a sample id')
parser$add_argument('--norm-method', type='character', default='NA',
    choices=c('RLE', 'CPM', 'TMM', 'upperquartile','NA'))
parser$add_argument('-o', '--output-file', type='character', required=TRUE,
    help='output file')
args <- parser$parse_args()

message('read count matrix: ', args$matrix)
mat <- read.table(args$matrix, header = TRUE, row.names=1, check.names=FALSE, sep='\t')

message('Read sample ids')
samples <- read.delim(args$sample_ids,sep="\n",stringsAsFactors=F,header=F)[,1]
message('Number of samples: ', length(samples))

mat <- as.matrix(mat[,samples])

norm_method <- args$norm_method

suppressPackageStartupMessages(library(edgeR))

if( norm_method == 'NA' ){
	message('Matrix output without normalization.')
	matrix_cpm <- mat
}else{
	message('Matrix output normalized by:',norm_method)
	matrix_cpm <- cpm(mat, method=norm_method)
}
message('Write subset and matrix to: ',args$output_file)
write.table(mat, args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
