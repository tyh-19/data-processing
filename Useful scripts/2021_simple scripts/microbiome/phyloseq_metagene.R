library(phyloseq)
setwd("C:/Users/Tao/Desktop/DESeq2 and phyloseq")

#otu table
otumat <- read.csv("counts-G_percent_forplot.csv",header = T)
head(otumat)
str(otumat)
#otumat <- as.matrix(otumat),2:1696 is genera number,2:7 is cacnertype
class(otumat)
otutable <- otumat[1:1695,2:7]
rownames(otutable) <- otumat[1:1695,1]
head(otutable)
str(otutable)
class(otutable)
otutable[is.na(otutable)] <- 0
otutable_mat <- as.matrix(otutable)
otutable_mat <- t(otutable_mat)
otutable_mat=apply(otutable_mat,1,as.numeric)
otutable_mat=apply(otutable_mat,2,as.numeric)

class(otutable_mat[1000,6])
rownames(otutable_mat) <- otumat[1:1695,1]
head(otutable_mat)

OTU = otu_table(otutable_mat, taxa_are_rows = TRUE)
OTU[is.na(OTU)] <- 0
OTU

#tax, 11865=7*1695
taxmat = matrix(sample(letters, 11865, replace = TRUE), nrow = nrow(otutable_mat), ncol = 7)
rownames(taxmat) <- rownames(otutable_mat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat
write.csv(taxmat,"temp.csv")
taxmat <- read.csv("taxmat.csv",header = T)

rownames(taxmat) <- rownames(otutable_mat)
taxmat <- taxmat[,-1]

TAX = tax_table(as.matrix(taxmat))
TAX

head(taxmat)
nrow(otutable_mat)
#phyloseq浣撶�?
physeq = phyloseq(OTU, TAX)
physeq

plot <- plot_bar(physeq, fill = "Phylum")
plot+xlab("")+ylab("Relative Abundance")+geom_bar(stat = 'identity',position = 'fill')+
  theme(axis.text=element_text(size=14,colour="black"),
        axis.text.x=element_text(angle=0,hjust = 0.5),
        axis.title =element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))
  
  #scale_fill_brewer(palette ="Set3",direction = 1)+
  

###heatmap

#sample data
sample_conc <- read.csv("C:/Users/Tao/Desktop/DESeq2 and phyloseq/sampledata_forplot.csv",header = TRUE)
head(sample_conc)
rownames(sample_conc) <- sample_conc[1:6,1]
class(sample_conc[4,2])
sampledata = sample_data(sample_conc)
class(sampledata)


#otu table
otumat <- read.csv("counts-G_cancertype_normalizeToMAX.csv",header = T)
head(otumat)
str(otumat)
#otumat <- as.matrix(otumat),2:1696 is genera number,2:7 is cacnertype
class(otumat)
otutable <- otumat[1:1695,2:7]
rownames(otutable) <- otumat[1:1695,1]
head(otutable)
str(otutable)
class(otutable)
otutable[is.na(otutable)] <- 0
otutable_mat <- as.matrix(otutable)
otutable_mat <- t(otutable_mat)
otutable_mat=apply(otutable_mat,1,as.numeric)
otutable_mat=apply(otutable_mat,2,as.numeric)

class(otutable_mat[1000,6])
rownames(otutable_mat) <- otumat[1:1695,1]
head(otutable_mat)

OTU = otu_table(otutable_mat, taxa_are_rows = TRUE)
OTU[is.na(OTU)] <- 0
OTU

#tax, 11865=7*1695
taxmat = matrix(sample(letters, 11865, replace = TRUE), nrow = nrow(otutable_mat), ncol = 7)
rownames(taxmat) <- rownames(otutable_mat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat
write.csv(taxmat,"temp.csv")
taxmat <- read.csv("taxmat_forplot.csv",header = T)

rownames(taxmat) <- rownames(otutable_mat)
taxmat <- taxmat[,-1]

TAX = tax_table(as.matrix(taxmat))
TAX

head(taxmat)
nrow(otutable_mat)
#phyloseq浣撶�?
physeq = phyloseq(OTU, TAX, sampledata)
physeq

data("GlobalPatterns")
gpac <- subset_taxa(physeq,Domain=="Virus")
gpac <- prune_taxa(names(sort(taxa_sums(gpac),TRUE)[1:5]), gpac)
plot_heatmap(gpac,low="#000033", high="#FF3300")

