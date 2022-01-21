##make phyloseq object
library(decontam)
library(phyloseq)

#sample data
sample_conc <- read.csv("C:/Users/Tao/Desktop/Lulab/2020/02.Project/02. Data/metagene/conc. of each sample.csv",header = TRUE)
head(sample_conc)
rownames(sample_conc) <- sample_conc[1:373,1]
class(sample_conc[4,2])
sampledata = sample_data(sample_conc)
class(sampledata)

#otu table
otumat <- read.csv("C:/Users/Tao/Desktop/Lulab/2020/02.Project/02. Data/metagene/counts-G.csv",header = FALSE)
head(otumat)
str(otumat)
#otumat <- as.matrix(otumat),2:1696 is genera number,2:374 is sample number
class(otumat)
otutable <- otumat[2:1696,2:374]
rownames(otutable) <- otumat[2:1696,1]
colnames(otutable) <- sample_conc[1:373,1]
head(otutable)
str(otutable)
class(otutable)
otutable[is.na(otutable)] <- 0
otutable_mat <- as.matrix(otutable)
otutable_mat <- t(otutable_mat)
otutable_mat=apply(otutable_mat,1,as.numeric)
otutable_mat=apply(otutable_mat,2,as.numeric)

class(otutable_mat[1000,372])
rownames(otutable_mat) <- otumat[2:1696,1]
head(otutable_mat)

OTU = otu_table(otutable_mat, taxa_are_rows = TRUE)
OTU[is.na(OTU)] <- 0
OTU

#tax, 11865=7*1695
taxmat = matrix(sample(letters, 11865, replace = TRUE), nrow = nrow(otutable_mat), ncol = 7)
rownames(taxmat) <- rownames(otutable_mat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat
setwd("C:/Users/Tao/Desktop")
write.csv(taxmat,"temp.csv")
taxmat <- read.csv("taxmat.csv",header = T)
rownames(taxmat) <- rownames(otutable_mat)
taxmat <- taxmat[,-1]

TAX = tax_table(as.matrix(taxmat))
TAX

nrow(taxmat)
nrow(otutable_mat)
#phyloseq体系
physeq = phyloseq(OTU, TAX, sampledata)
physeq

##Identify Contaminants - Frequency
contamdf.freq <- isContaminant(physeq, method="frequency", conc="Library_yield.ng.")
head(contamdf.freq)
write.csv(contamdf.freq,"C:/Users/Tao/Desktop/Lulab/2020/02.Project/02. Data/metagene/decontam_G.csv")
