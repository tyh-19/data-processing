# first input .rnk and .gmt as rnk.file, gmx.file
rnk.file <- "test.rnk"
gmt.file <- "test.gmx"

# extract ranklist and geneset
ranks <- read.table(rnk.file, header=T, colClasses = c("character", "numeric"))

ranks <- setNames(ranks$Rank,ranks$GeneID)

pathways <- read.table(gmt.file)
pathways <- t(pathways)
write.table(pathways,file="test.gmt",col.names = F,row.names = F, quote = F, sep = "\t")
pathways <- gmtPathways("test.gmt")

# check your ranklist and geneset
str(ranks)
str(head(pathways))

# run user-defined GSEA. 
# Question about warning message, refer to https://www.gitmemory.com/issue/ctlab/fgsea/79/703686364
fgseaRes <- fgsea(pathways, ranks, minSize=0, maxSize=500, nperm = 1000)
write.csv(fgseaRes,"fgsea.csv")

head(fgseaRes[order(pval), ])

files <- as.data.frame(dir())
files2 <- files[grep("CRC-",files$`dir()`, value = F),]
files2 <- as.data.frame(files2)
files2[2,]
library(ggpubr)
files2[1,]
i <- 1
while (i+1 < nrow(files)){
gfoldcount <- read.table(as.character(files2[i,]))

featurecount <- read.table(as.character(files2[i+2,]))

f <- data.frame(featurecount$V1,featurecount$V3)
f2 <- f[order(f[,1]),]

g <- data.frame(gfoldcount$V1,gfoldcount$V3)

cor <- data.frame(g$gfoldcount.V3,f2$featurecount.V3)

head(cor)

p <- ggplot(data=cor, aes(x=g.gfoldcount.V3,y=f2.featurecount.V3))+geom_point(color="red")+
  stat_cor(data=cor, method = "pearson")+
  xlab("gfoldcount")+
  ylab("featurecount")

ggsave(p,filename = paste0(as.character(files2[i,]),".pdf"))

i <- i+2
}



