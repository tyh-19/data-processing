# **Data processing pipeline**

### 01. Raw data

> paired end sequencing data 

### 02. Quality control: fastqc

```sh
#!/bin/bash
input=$"/BioII/lulab_b/wangsiqi/exRNA/exRNA-panel/pico/01.rawdata/"
output=$"/home/taoyuhuan/liquid/STAD/fastqc"

files=$(ls ${input}|grep STAD|grep 1.fq.gz|cut -d "_" -f 1)
for filename in ${files} 
do
	mkdir ${output}/${filename}
	echo ${filename} " is processing,please wait..."
     	fastqc ${input}/${filename}_1.fq.gz ${input}/${filename}_2.fq.gz --outdir ${output}/${filename} --noextract       	
done
```

### 03. Trim: cutadaptor & trimGC.py

+ cutadaptor

```sh
#!/bin/bash
input="/home/taoyuhuan/liquid/STAD/01.raw"
output="/home/taoyuhuan/liquid/STAD/03.trim"
adaptor1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adaptor2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

files=$(ls ${input} | grep _1.fq.gz | cut -d "_" -f 1,2)

for filenames in ${files}
	do
	echo ${filenames}
	cutadapt --pair-filter any --trim-n --match-read-wildcards --minimum-length 15 -q 10,10 -a ${adaptor1} -A ${adaptor2} -o ${output}/${filenames}_1.fq.gz -p ${output}/${filenames}_2.fq.gz ${input}/${filenames}_1.fq.gz ${input}/${filenames}_2.fq.gz > ${output}/${filenames}_cutAdaptor.log
done
```

+ trimGC.py

> in order to trim several template switch step added GGG(CCC)

```sh
#!/bin/bash
#require trimGC.py [write by SBB, revised by JYF,CYH; please put them in the same dictionary or use full path]
input="/home/taoyuhuan/liquid/STAD/03.trim"
output="/home/taoyuhuan/liquid/STAD/04.TSO"

file=$(ls ${input} | grep _1.fq.gz | cut -d "_" -f 1,2)

for filenames in ${file}
	do
	echo ${filenames}
	python trimGC.py -i ${input}/${filenames} -o ${output}/${filenames} -m 15 -s reverse > ${output}/${filenames}_TSO.log
done
```

### 04. Mapping: [STAR](https://www.bioinfo-scrounger.com/archives/288/)

```sh
#!/bin/bash
input="/home/taoyuhuan/liquid/STAD/04.TSO"
output="/home/taoyuhuan/liquid/STAD/05.mapping"

files=$(ls ${input}|grep 1.fq.gz|cut -d "_" -f 1,2)
for filename in ${files} 
do
	cd ${output}
	echo ${filename}
       	STAR --runThreadN 6 \
	     --runMode alignReads \
	     --readFilesCommand gunzip \
	     --outFileNamePrefix ${filename} \
         --genomeDir /BioII/lulab_b/shared/genomes/human_hg38/index/STAR_hg38_index \
         --readFilesIn ${input}/${filename}_1.fq.gz ${input}/${filename}_2.fq.gz \
         --outSAMtype BAM SortedByCoordinate        	
done
```

###  05. Featurecount:  FeatureCounts

```sh
#!/bin/bash
input="/home/taoyuhuan/liquid/STAD/05.mapping_fq"
output="/home/taoyuhuan/liquid/STAD/06.featurecount"

for file in $(ls ${input}|grep Aligned.sortedByCoord.out.bam)
	do
	echo ${file}
	featureCounts -T 6 -s 2 -p -t exon -g gene_id -a /home/taoyuhuan/reference/hg38_annotation.gtf -o ${output}/${file}.featurecount.txt ${input}/${file}
done
```

### 06. Build matrix

```sh
#! /bin/bash
input="/home/taoyuhuan/liquid/STAD/06.featurecount"
output="/home/taoyuhuan/liquid/STAD/07.matrix"

first=$(ls ${input}|grep .txt|grep -v .summary|head -n 1)
file=$(ls ${input}|grep .txt|grep -v .summary|cut -d "." -f 1) 
echo ${first}
#get gene id, delete first line(featurecount step generated)
cut -f 1 ${input}/${first} > ${output}/firstcol.tmp.txt
sed -i '1d' ${output}/firstcol.tmp.txt
#get featurecount of each group, delete first line and alter second line to group info
for filename in ${file}
	do
	echo ${filename}
	cut -f 7 ${input}/${filename}.sortedByCoord.out.bam.featurecount.txt > ${output}/${filename}.tmp.txt
	sed -i '1,2d' ${output}/${filename}.tmp.txt
	#get the group info from filename. Since  cut only process file or dir, i touch a temp.txt to process filename.
	touch ${filename}.txt
	group=$(ls ${filename}.txt|cut -d "_" -f 1)
	echo ${group}
	sed -i "1i \\${group}" ${output}/${filename}.tmp.txt
	rm ${filename}.txt
done
#merge to a matrix
cd ${output}
paste *.tmp.txt > Matrix_rawcount.txt

rm ${output}/*.tmp.txt
```

### 07. Different gene expression in R

+ DESeq

```R
##载入DESeq2
> library(DESeq2)

##读入matrix的行名（transcriptID）
> rowname <-read.csv("/home/taoyuhuan/liquid/STAD/07.matrix/rowname.csv",header=TRUE)

##读入matrix内的rawcounts
> rawcounts <-read.csv("/home/taoyuhuan/liquid/STAD/07.matrix/Matrix_rawcount.csv",header=TRUE)

##将rowname转化为int型，并作为rawcounts的rowname，将新的有transcriptID作行名的matrix输出保存
> rowname_int <- gsub("\\.\\d*", "",rowname[,1])
> row.names(rawcounts) <- rowname_int
> write.csv(x=rawcounts,file="/home/taoyuhuan/liquid/STAD/08.DEseq/rawcounts.csv")

##制作列处理矩阵（coldata）
> condition <- factor(c("NC","NC","STAD","STAD"))
> coldata <- data.frame(row.names=colnames(rawcounts),condition)
> coldata
       condition
NC            NC
NC.1          NC
STAD        STAD
STAD.1      STAD

##从matrix制作DESeqDataSet，标准化（命令：DESeq（dds））后，查看
> dds <- DESeqDataSetFromMatrix(rawcounts, DataFrame(condition), design= ~ condition )
> head(dds)
class: DESeqDataSet
dim: 6 4
metadata(1): version
assays(1): counts
rownames(6): ENST00000456328 ENST00000450305 ... ENST00000473358
  ENST00000469289
rowData names(0):
colnames(4): NC NC.1 STAD STAD.1
colData names(1): condition
> dds_normalized <- DESeq(dds)
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
> resultsNames(dds_normalized)
[1] "Intercept"            "condition_STAD_vs_NC"

##将标准化后的结果输出并命名为res
> res <- results(dds_normalized)
> head(res)
log2 fold change (MLE): condition STAD vs NC
Wald test p-value: condition STAD vs NC
DataFrame with 6 rows and 6 columns
                         baseMean     log2FoldChange            lfcSE
                        <numeric>          <numeric>        <numeric>
ENST00000456328  27.1154375587156  -1.81571222116039 2.13598414043766
                              stat            pvalue              padj
                         <numeric>         <numeric>         <numeric>
ENST00000456328 -0.850058849588815 0.395292368471226 0.461309011695085
> mcols(res,use.names= TRUE) # 查看res矩阵每一列的含义
DataFrame with 6 rows and 2 columns
                       type                                  description
                <character>                                  <character>
baseMean       intermediate    mean of normalized counts for all samples
log2FoldChange      results log2 fold change (MLE): condition STAD vs NC
lfcSE               results         standard error: condition STAD vs NC
stat                results         Wald statistic: condition STAD vs NC
pvalue              results      Wald test p-value: condition STAD vs NC
padj                results                         BH adjusted p-values
> summary(res)
out of 27955 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 3849, 14%
LFC < 0 (down)     : 1183, 4.2%
outliers [1]       : 0, 0%
low counts [2]     : 15586, 56%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
#上调基因3849，下调基因1183个
> table(res$padj<0.05) # 取padj小于0.05的数据，得到3709行
FALSE  TRUE
 8660  3709

##将标准化后的counts和差异表达结果合并成一个dataframe，按照padj从小到大排序，取出padj<0.5，变化倍数>2的transcriptID，并取出rawcounts（前4列，根据绘制的rawcount matrix），输出保存
> cbind <- cbind(as.data.frame(counts(dds_normalized,normalize=TRUE)),results(dds_normalized))
> head(cbind)
DataFrame with 6 rows and 10 columns
                               NC             NC.1            STAD
                        <numeric>        <numeric>       <numeric>
ENST00000456328  27.5769208422463  56.600979086782 21.777261664733
                          STAD.1          baseMean     log2FoldChange
                       <numeric>         <numeric>          <numeric>
ENST00000456328 2.50658864110104  27.1154375587156  -1.81571222116039
                           lfcSE               stat            pvalue
                       <numeric>          <numeric>         <numeric>
ENST00000456328 2.13598414043766 -0.850058849588815 0.395292368471226
                             padj
                        <numeric>
ENST00000456328 0.461309011695085

> cbind_ordered <- cbind[order(cbind$padj),]
> cbind_ordered_diffgene <- subset(cbind_ordered,padj<0.5,log2FoldChange > 1 | log2FoldChange < -1)
> write.csv(cbind_ordered_diffgene,"/home/taoyuhuan/liquid/STAD/08.DEseq/cbind_ordered_diffgene_padj0.5_foldchange2.csv")
> diffgene_STAD2vsNC2 <- cbind_ordered_diffgene[,1:4]
> head(diffgene_STAD2vsNC2)
DataFrame with 6 rows and 4 columns
                              NC             NC.1             STAD
                       <numeric>        <numeric>        <numeric>
ENST00000643151 2.12130160324971   4.353921468214 172294.435204146
                          STAD.1
                       <numeric>
ENST00000643151 154027.365407018
> write.csv(diffgene_STAD2vsNC2,"/home/taoyuhuan/liquid/STAD/08.DEseq/diffgene_STAD2vsNC2.csv")
```

+ heatmap

```R
> require(pheatmap)
Loading required package: pheatmap
> pdf("/home/taoyuhuan/liquid/STAD/08.DEseq/pheatmap_STAD2vsNC2.pdf",pointsize=1)
> dev.list()
pdf
  2
> pheatmap(diffgene_STAD2vsNC2)
> dev.off()

> diffgene_STAD2vsNC2_30 <- diffgene_STAD2vsNC2[1:30,]
> write.csv(diffgene_STAD2vsNC2_30,"/home/taoyuhuan/liquid/STAD/08.DEseq/diffgene_STAD2vsNC2_30.csv")
> pdf("/home/taoyuhuan/liquid/STAD/08.DEseq/pheatmap_STAD2vsNC2.pdf",pointsize=1)
> pheatmap(diffgene_STAD2vsNC2_30, scale = "row")
> dev.off()
```

<div align=left><img src="C:\Users\Tao\Desktop\data processing\heatmap.png" style="zoom:40%;" /><div>

+ TranscriptID > name change by biomaRt in R

```R
> library("biomaRt")
##仅输入需要转换的list，不要输入header
> temp <- read.csv("/home/taoyuhuan/liquid/STAD/08.DEseq/CONVERT.csv",header=F)
> head(temp)
##导入一个识别的库
> mart <- useMart("ensembl","hsapiens_gene_ensembl") 
##仅将transcriptID转换为external_gene_name,并同时保留
> temp1 <- getBM(attributes=c('ensembl_transcript_id','external_gene_name'), filters= 'ensembl_transcript_id', values = temp, mart = mart)
> head(temp1)
##看看可以转化成其他的什么，然后以transcriptID为filter，进行多类型的同时注释，例如gene_source/description，并保存
head(listAttributes(mart))
temp1 <- getBM(attributes=c('ensembl_transcript_id','external_gene_name','external_gene_source','description'), filters= 'ensembl_transcript_id', values = temp, mart = mart)
head(temp1)
write.csv(temp1,"converted.csv")
```

<div align=left><img src="C:\Users\Tao\Desktop\data processing\transcriptID to annotation.png" style="zoom:40%;" /><div>

> 在PPT中将heatmap 和 gene name list拼接



### 08. GO & KEGG analysis 

+ GO and bar plot

> 取出transcriptID（ensembl），至DAVID 网页分析工具分析，得到MF/BP/CC三类注释，作bar plot

```R
> library(ggplot2)
> GO_BP <- read.table("/home/taoyuhuan/liquid/STAD/08.DESeq/diffgene_GO.txt",header = T,sep="\t")
> pdf("/home/taoyuhuan/liquid/STAD/08.DESeq/GO_BP_2.pdf", pointsize = 1)
> ggplot(data=GO_BP)+
  geom_bar(aes(x=reorder(Term,Count),y=Count, fill=-log10(PValue)), stat='identity') + 
  coord_flip() +
  scale_fill_gradient(expression(-log["10"](P.value)),low="red", high = "blue") +
  xlab("") +
  ylab("Gene count")
> dev.off()
```

<div align=left><img src="C:\Users\Tao\Desktop\data processing\GO_BP.png" style="zoom:40%;" /><div>



+ KEGG and bubble plot

> 取出transcriptID（ensembl），至DAVID 网页分析工具分析，得到KEGG结果，作bubble plot

```R
> library(ggplot2)
> KEGG <- read.csv("diffgene_NC2vsSTAD2_KEGG_top10.csv",header=T)

> pdf("KEGG2.pdf",pointsize=1)
> ggplot(KEGG,aes(x=Fold.Enrichment,y=Term))+
+ geom_point(aes(size=Count,color=-1*log10(PValue)))+
+   scale_colour_gradient(low="green",high="red")+
+   labs(
+        color=expression(-log[10](P.value)),
+        size="Gene number",
+        x="Fold enrichment"
+        # y="Pathway name",
+        # title="Pathway enrichment")
+       )+
+   theme_bw()+
+   theme(
+     axis.text.y = element_text(size = rel(1.3)),
+     axis.title.x = element_text(size=rel(1.3)),
+     axis.title.y = element_blank()
+   )
Warning message:
Removed 9 rows containing missing values (geom_point).
> dev.off()
```

<div align=left><img src="C:\Users\Tao\Desktop\data processing\KEGG.png" style="zoom:40%;" /><div>

