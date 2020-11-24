Copy commands may not work due to format change between Mac-linux, Win-linux, type it. 

# General

## txt manipulate

### colum joint and count col

```bash
# we can set delimited by -d ''
paste file1 file2
# count first line words
head -1 file | wc -w
```

### Convert dos to unix (remove ^M)

```bash
# -n means write plasma_feature.txt to plasma in unix mode
dos2unix -n plasma_featurecount.txt plasma.txt
# -o means write to old file
dos2unix -o plasma_featurecount.txt
```

### delete first col

```bash
cut -d "," -f2- featurecount_summary.csv > PBMC_rawcount.csv
```

### sum col

```bash
cat CRC-2418658-T-circRNA | awk '{sum+=$2} END{print "Sum = ", sum}'
```

### change element by tr

```bash
# tr 'x' 'y' < file
cat /BioII/lulab_b/taoyuhuan/liquid/GSEA/expression-counts.txt | head -n 1 | tr '\t' '\n' > /BioII/lulab_b/taoyuhuan/liquid/GSEA/all.txt
```

### Windows to linux

```bash
## windows will add ^M to the end of line before $, while linux only have $ at the end
cat -A file | while read line
do
id=`echo $line | cut -d"^" -f1`
# id is the same end with linux
```

### Intersect, union and except

```bash
## intersect
cat a.txt b.txt | sort | uniq -d > intersected.txt
## union
cat a.txt b.txt | sort | uniq > union.txt
##
cat a.txt b.txt | sort | uniq -u > except.txt
```

### sed

```bash
# delete line contian special string, change raw file (-i), default do not change raw
sed -i '/baidu.com/d' domain.file
# delete special string in each site, retain other
sed 's/_Abundance//g' Pathway_Abundance_summary.tsv > test.tsv
# delete null line
sed -i '/^$/d' a.txt
```

### grep

```
# show other lines without abc
grep -v abc
```

### awk&sort

```bash
### filter colum2 != 0, sort by colum2(-k 2,sep by '\t'), from high to low(-r) as number(-n)
cat file | awk -F '\t' '$2 != 0.0' | sort -n -r -k 2
```

### Judge Phred+33/64(not suit for sra transformed fastq)

```bash
zcat /BioII/lulab_b/jinyunfan/projects/exSEEK/exSeek-dev/output/pico_reverse/unmapped/CRC-2124325/circRNA_1.fastq.gz | head -100 | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding!";}'
```

### A. SeqPrep

+ Installation
  + Download [SeqPrep](https://codeload.github.com/jstjohn/SeqPrep/zip/master)
  + cd /path/to/SeqPrep
  + make install

+ Usage

  ```bash
  /home/taoyuhuan/tools/SeqPrep-master/SeqPrep
  
  /home/taoyuhuan/tools/SeqPrep-master/SeqPrep -f circRNA_1.fastq.gz -r circRNA_2.fastq.gz -1 CRC-2124325_1.fastq.gz -2 CRC-2124325_2.fastq.gz -3 CRC-2124325_1_discard.fastq.gz -4 CRC-2124325_2_discard.fastq.gz -s CRC-2124325.fastq.gz
  ```
  
+ 

# RNA different expression

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

###  05. Featurecount:  [FeatureCounts](http://bioinf.wehi.edu.au/featureCounts/)

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

+ [DESeq](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

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

<div align=left><img src="C:\Users\Tao\Documents\GitHub\data processing\heatmap.png" style="zoom:40%;" /><div>


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

<div align=left><img src="C:\Users\Tao\Documents\GitHub\data processing\transcriptID to annotation.png" style="zoom:40%;" /><div>


> 在PPT中将heatmap 和 gene name list拼接



### 08. Sample Correlation

```R
> install.packages("pheatmap")
> stad_raw=read.csv("C:/Users/Tao/Desktop/rawcounts.csv",header=T,row.names=1)
> head(stad_raw)
                 NC NC.1 STAD STAD.1
ENST00000456328 182  104    6      1
ENST00000450305   0    0    0      0
ENST00000488147   2   18    0      2
ENST00000619216   0    2    0      0
ENST00000473358   0    0    0      0
ENST00000469289   0    0    0      0
> cor_stad_raw_matrix=cor(stad_raw)
> head(cor_stad_raw_matrix)
              NC      NC.1      STAD    STAD.1
NC     1.0000000 0.9792599 0.2784779 0.2655405
NC.1   0.9792599 1.0000000 0.2169506 0.2109196
STAD   0.2784779 0.2169506 1.0000000 0.9769124
STAD.1 0.2655405 0.2109196 0.9769124 1.0000000
> write.csv(cor_stad_raw_matrix,'cor_stad_raw_matrix.csv')
> pheatmap(cor_stad_raw_matrix,cluster_rows=F,cluster_cols=F,display_numbers=T)
```



 <img src="C:\Users\Tao\Documents\GitHub\data processing\cor_raw.png" alt="cor_rawcounts" style="zoom:70%;" /><img src="C:\Users\Tao\Documents\GitHub\data processing\cor_normalizedbyDESeq.png" alt="cor_normalizedbyDEseq" style="zoom:70%;" />

### 09. Volcano plot

```R
diffgene <- read.table("/home/test/share/RNAseq_homework/gene_exp_homework.diff",header=T,sep="\t")


diffgene$threshold <- as.factor(ifelse(diffgene$q_value < 0.05 & abs(diffgene$log2.fold_change.) >=1,ifelse(diffgene$log2.fold_change. > 1 ,'Up','Down'),'Not'))
volcanoplot_gene <- ggplot(data=diffgene, aes(x=log2.fold_change., y =-log10(q_value), color=threshold,fill=threshold,label=gene)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_point(size=1) +
  xlim(c(-3, 3)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
  theme(legend.position="right",
  panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=12),
        axis.text.y = element_text(face="bold",  color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(x="log2(fold_change)",y="-log10 (q_value)",title="Volcano plot of different expression gene(geneID)", face="bold")

volcanoplot_label_gene <- volcanoplot_gene + geom_text_repel(force=0.01,nudge_y=0.05) + labs(title = "Volcanoe plot of DEG(gene)")

ggsave("/home/test/share/RNAseq_homework/volcano_plot_tyh_genelabe_repel.pdf", plot=volcanoplot_gene, height = 10, width = 10)
```



### 10. GO & KEGG analysis 

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

<div align=left><img src="C:\Users\Tao\Documents\GitHub\data processing\GO_BP.png" style="zoom:40%;" /><div>




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

<div align=left><img src="C:\Users\Tao\Documents\GitHub\data processing\KEGG.png" style="zoom:40%;" /><div>
### 11. Gfold

```bash
# convert featurecount results to gfold read_cnt format
# colnames of read_cnt: Gene Symbol, Gene_name, Readcount, Gene_exon_length, RPKM
# we can fake Gene_name and RPKM while not influence results
cat /BioII/lulab_b/taoyuhuan/liquid/CRC_piared_tissue/09.featureCounts/CRC-2211756-N | awk '{print $1,"GeneName "$7,$6,"1"}' | tr ' ' '\t' > test.txt; sed -i '1,2d' test.txt
```

### 10.Quality_Control

Get mapping ratio in Log.final.out

```bash
grep Uniquely mapped reads CRC-2211756-N-NOspikein-NOunivec-NOrRNA/*Log.final.out
```

Get reads number in bam (remove duplicated)

```bash
# !/bin/bash

input="/BioII/lulab_b/jinyunfan/projects/exSEEK/exSeek-dev/output/shaozhen-1007/bam/"
output="/BioII/lulab_b/taoyuhuan/liquid/Quality_control/summary"
target="UniVec"

file=$(ls ${input} | cut -d "-" -f 1,2,3)

echo 'Sample ID' > ${output}/${target}_ID.txt
echo ${target} > ${output}/${target}.txt

for filename in ${file}
do
echo ${filename} >> ${output}/${target}_ID.txt
samtools view -c -f 0x2 ${input}/${filename}/${target}.bam >> ${output}/${target}.txt
done

paste ${target}_ID.txt ${output}/${target}.txt > ${output}/${target}_readcount.txt
```

### 11.Gfold

Convert featurecount result to gfold format

```bash
cat /BioII/lulab_b/taoyuhuan/liquid/CRC_piared_tissue/09.featureCounts/CRC-2211756-N | awk '{print $1,"GeneName "$7,$6,"1"}' | tr ' ' '\t' > test.txt; sed -i '1,2d' test.txt
```

Gfold count

```bash
#!/bin/bash

input="/BioII/lulab_b/taoyuhuan/liquid/CRC_piared_tissue/03.4.mappinghg38_cnode/Bam_hg38"
output="/BioII/lulab_b/taoyuhuan/liquid/CRC_piared_tissue/10.DE/1vs1_gfold"
ref="/BioII/lulab_b/shared/genomes/hg38/gtf/long_RNA.gtf"

file=$(ls ${input} | grep .bam | grep '2416785\|2418277\|2418488\|2418503\|2418658' | grep -v plasma | cut -d "-" -f 1,2,3)
for filename in ${file}
do
echo ${filename}
samtools view ${input}/${filename}-NOspikein-NOunivec-NOrRNA_Aligned_hg38.bam | gfold count -ann ${ref} -tag stdin -o ${output}/${filename}.read_cnt
done
```

Gfold diff

```bash
# input featurecount and gfold count respectively, then we can compare the difference between two count method
#!/bin/bash

input="/BioII/lulab_b/taoyuhuan/liquid/CRC_piared_tissue/10.DE/1vs1_gfold/featurecount2gfoldcount"
output="/BioII/lulab_b/taoyuhuan/liquid/CRC_piared_tissue/10.DE/1vs1_gfold/diff"


file=$(ls ${input} | grep '2416785\|2418277\|2418488\|2418503\|2418658' | grep T | cut -d "." -f 1 | cut -d "-" -f 1,2)
for filename in ${file}
do
echo ${filename}
gfold diff -acc T -s1 ${input}/${filename}-T -s2 ${input}/${filename}-N -suf .txt -o ${output}/${filename}-T_vs_N_featurecount.diff
done
```

Gfold for metabolism

```bash

```



### 12.TCGA

Countsummary

```bash
#!/bin/bash
# COADcount path=/readonly/Share2/home/lulab1/TCGA/processed/Gene_Expression_Quantification/RNA-seq/COAD/htseq_counts
# READcount path=/readonly/Share2/home/lulab1/TCGA/processed/Gene_Expression_Quantification/RNA-seq/READ/htseq_counts

input=$1

colnames=$(ls $1 | grep .gz | cut -d "." -f 1)

echo ${colnames}
mkdir temp
rm /BioII/lulab_b/taoyuhuan/liquid/GSEA/TCGA_CRC_698/temp/*_temp.txt

rowname=$(ls $1 | grep .gz | head -1)
echo "ID        ${colname}" > ./temp/0000000000000000rowname_temp.txt
zcat $1/${rowname} | awk -F '\t' '{print $1}' >> ./temp/0000000000000000rowname_temp.txt

for colname in ${colnames}
do
echo "ID        ${colname}" > ./temp/${colname}.txt
zcat $1/${colname}.htseq.counts.gz >> ./temp/${colname}.txt
cut -d "        " -f2- ./temp/${colname}.txt > ./temp/${colname}_temp.txt
done

paste ./temp/*_temp.txt > htseqcount_summary_READ.txt
```

Group info

```bash
# check all groups
awk -F '   ' '{print $29}' COAD_RNA_Seq_summary.tsv | SORT -u | uniq

# get specific ID in each group
cat COAD_RNA_Seq_summary.tsv | grep "Recurrent Solid Tumor" | awk -F '     ' '{print $5}' | grep htseq | cut -d '.' -f 1 > COAD_RecurrentSolidTumor.txt

# check if all ID grouped
(base) [taoyuhuan@cnode metadata]$ wc -l *.txt
    1 COAD_Metastatic.txt
  478 COAD_PrimarySolidTumor.txt
    1 COAD_RecurrentSolidTumor.txt
   41 COAD_SolidTissueNormal.txt
  521 total
(base) [taoyuhuan@cnode metadata]$ cat COAD_RNA_Seq_summary.tsv | grep htseq.counts.gz | wc -l
	521
```

### 13.bigwig normalize

```bash

```

### 14.GSEA

+ Make genesets by DE

  ```bash
  # get inhoused CRC featurecount summarized matrix and TCGA READ/COAD matrix
  
  # get group information
  cat featurecount_summary.txt | head -1 | tr '    ' '\n' | grep plasma | grep -v "2399129\|2384058" > CRC_plasma.txt
  cat featurecount_summary.txt | head -1 | tr '    ' '\n' | grep NC > NC_plasma.txt
  cat featurecount_summary.txt | head -1 | tr '    ' '\n' | grep NC | grep PKU > NC_PKU_plasma.txt
  cat featurecount_summary.txt | head -1 | tr '    ' '\n' | grep CRC | grep T > CRC_Tumor.txt
  cat featurecount_summary.txt | head -1 | tr '    ' '\n' | grep CRC | grep N > CRC_Normal.txt
  
  # DE (with replicates)
  Rscript ../scripts/differential_expression.R -i ./raw_matrix/featurecount_summary.txt -p ./group/CRC_Tumor.txt -n ./group/CRC_Normal.txt -m edger_glmlrt -o CRC_Tumor_vs_Normal.txt
  
  Rscript ../scripts/differential_expression.R -i ./raw_matrix/htseqcount_summary_COAD.txt -p ./group/COAD_PrimarySolidTumor.txt -n ./group/COAD_SolidTissueNormal.txt -m edger_glmlrt -o ./DE/TCGA_Tumor_vs_Normal.txt
  
  Rscript ../scripts/differential_expression.R -i ./raw_matrix/htseqcount_summary_READ.txt -p ./group/READ_PrimarySolidTumor.txt -n ./group/READ_SolidTissueNormal.txt -m edger_glmlrt -o ./DE/READ_Tumor_vs_Normal.txt
  
  # DE (without replicates)
  # Convert featurecount2gfold format
  # !/bin/bash
  input="/BioII/lulab_b/taoyuhuan/liquid/CRC_piared_tissue/09.featureCounts_repair"
  output="/BioII/lulab_b/taoyuhuan/liquid/GSEA/Geneset_1027/gfold"
  file=$(ls ${input} | grep -v .summary )
  for filename in ${file}
  do
  echo ${filename}
  cat ${input}/${filename} | awk '{print $1,"GeneName "$7,$6,"1"}' | tr ' ' '\t' > ${output}/${filename}_gfold.txt; sed -i '1,2d' ${output}/${filename}_gfold.txt
  done
  
  # Gfold diff
  
  # input featurecount and gfold count respectively, then we can compare the difference between two count method
  # !/bin/bash
  
  input="/home/taoyuhuan/liquid/GSEA/Geneset_1027/gfold"
  output="/home/taoyuhuan/liquid/GSEA/Geneset_1027/gfold/diff_forGMX"
  
  
  file=$(ls ${input} | grep T | cut -d "." -f 1 | cut -d "-" -f 1,2)
  for filename in ${file}
  do
  echo ${filename}
  #gfold diff -acc T -s1 ${input}/${filename}-N_gfold -s2 ${input}/${filename}-T_gfold -suf .txt -o ${output}/${filename}-T_vs_N_featurecount.diff
  
  cat ${output}/${filename}-T_vs_N_featurecount.diff | awk '!($3==0){print $1,$3}' | tr ' ' '\t' > ${input}/geneset/${filename}_T_vs_N.tmp
  sed -i '1,9d' ${input}/geneset/${filename}_T_vs_N.tmp
  cat ${input}/geneset/${filename}_T_vs_N.tmp | tail -n +2 | sort -r -t ' ' -n -k 2 > ${input}/geneset/${filename}_T_vs_N.txt
  sed -i '1i\GeneID       Gfold' ${input}/geneset/${filename}_T_vs_N.txt
  rm ${input}/geneset/${filename}_T_vs_N.tmp
  done
  
  # merge to .gmx
  #!/bin/bash
  input_gfold="/home/taoyuhuan/liquid/GSEA/Geneset_1027/gfold/geneset"
  input_DE="/home/taoyuhuan/liquid/GSEA/Geneset_1027/DE"
  output="/home/taoyuhuan/liquid/GSEA/Geneset_1027/geneset"
  
  file=$(ls ${input_gfold} | grep .txt | grep "2416785\|2418277\|2418488\|2418503\|2418658" | cut -d '.' -f 1)
  
  # sort is from low to high, head is down while tail is up
  for filename in ${file}
  do
  echo ${filename}
  cat ${input_gfold}/${filename}.txt | tail -n +2 | sort -t "     " -n -k 2 | head -500 | awk '{print $1}'> ${output}/${filename}_Down_500.txt
  cat ${input_gfold}/${filename}.txt | tail -n +2 | sort -t "     " -n -k 2 | tail -500 | awk '{print $1}'> ${output}/${filename}_Up_500.txt
  sed -i "1i/${filename}_Up_500" ${output}/${filename}_Up_500.txt
  sed -i "1i/${filename}_Down_500" ${output}/${filename}_Down_500.txt
  done
  
  file2=$(ls ${input_DE} | grep .txt | cut -d '.' -f 1)
  
  for filename2 in ${file2}
  do
  echo ${filename2}
  cat ${input_DE}/${filename2}.txt | tail -n +2 | awk '($6<=0.05){print $0}' | sort -t "  " -n -k 2 | head -500 | awk '{print $1}'> ${output}/${filename2}_Down_500.tmp
  cat ${input_DE}/${filename2}.txt | tail -n +2 | awk '($6<=0.05){print $0}' | sort -t "  " -n -k 2 | tail -500 | awk '{print $1}'> ${output}/${filename2}_Up_500.tmp
  sed -i "1i/${filename2}_Up_500" ${output}/${filename2}_Up_500.tmp
  awk -F '|' '{print $1}' ${output}/${filename2}_Up_500.tmp > ${output}/${filename2}_Up_500.txt
  sed -i "1i/${filename2}_Down_500" ${output}/${filename2}_Down_500.tmp
  awk -F '|' '{print $1}' ${output}/${filename2}_Down_500.tmp > ${output}/${filename2}_Down_500.txt
  rm ${output}/*.tmp
  done
  
  paste ${output}/*500.txt > ${output}/Geneset_1027.gmx
  
  rm ${output}/*.txt
  
  ```

  

+ Make ranksets by ratio

  ```bash
  # for comparisom between group
  python ../scripts/TPM.py -i featurecount_summary.txt -o TPM.txt
  
  # filter genes TPM < 1, showup <50%
  python ../scripts/filter.py -i TPM.txt -o Filtered_TPM.txt -p 0.5 -s 1
  
  calculate average in python:https://jingyan.baidu.com/article/4ae03de3b436233eff9e6b94.html
  # calculate group average
  ## All NC average (N=99)
  python ../scripts/row_average.py -i Filtered_TPM.txt -g NC -u CRC -o group_average.txt
  ## All PKU NC average (N=54)
  python ../scripts/row_average.py -i group_average.txt -g NC_PKU -u CRC -o group_average.txt
  ## All ShH NC average (N = 33)
  python ../scripts/row_average.py -i group_average.txt -g NC_ShH -u CRC -o group_average.txt
  ## All ChQ NC average (N = 12)
  python ../scripts/row_average.py -i group_average.txt -g NC_ChQ -u CRC -o group_average.txt
  ## All CRC-plasma average (N = 84)
  python ../scripts/row_average.py -i group_average.txt -g CRC -u N,T,plasma,PBMC -o group_average.txt
  ## All paired Normal (N =15)
  python ../scripts/row_average.py -i group_average.txt -g N -u NC,plasma,PBMC,T -o group_average.txt
  ## All paired Tumor (N = 15)
  python ../scripts/row_average.py -i group_average.txt -g T -u plasma,PBMC,N -o group_average.txt
  ```

  + row_average.py

  ```python
  import warnings
  warnings.simplefilter(action='ignore', category=FutureWarning)
  import numpy as np
  import pandas as pd
  import argparse
  import sys
  parser = argparse.ArgumentParser(description='Calculate specific row average')
  parser.add_argument('--input', '-i', type=str, required=True, help='input count matrix')
  parser.add_argument('--output','-o',type=str,required=True,help='output matrix')
  parser.add_argument('--group','-g',type=str,required=True,help='interested group')
  parser.add_argument('--ungroup','-u',type=str,required=True,help='uninterested group, default remove "Average_". If more than one element, split by ,')
  args = parser.parse_args()
  
  def in_uninterested(n):
          i=0
          j=0
          while i < len(ungroup_list):
                  #print(ungroup_list[i])
                  if ungroup_list[i] in n:
                          i=i+1
                          j=j+1
                  else:
                          i=i+1
          if j == 0:
                  return 1
          else:
                  return 0
  
  print("Load data ...")
  df = pd.read_csv(args.input,index_col=0,sep="\t")
  
  ungroup_list = ['Average_']
  ungroup_list.extend(args.ungroup.split(','))
  
  print("Number of intersted:")
  print(len(list(df.filter(like=args.group).columns)))
  print((list(df.filter(like=args.group).columns)))
  
  # axis = 1, row elements average; while axis = 0, column elements average
  df['Average_'+args.group] = df[list(filter(in_uninterested,list(df.filter(like=args.group).columns)))].mean(axis=1)
  
  print("Unwanted list:")
  print(ungroup_list)
  print("Number of interested after remove unwanted:")
  print(len(list(filter(in_uninterested,list(df.filter(like=args.group).columns)))))
  print("Average of:")
  print((list(filter(in_uninterested,list(df.filter(like=args.group).columns)))))
  df.to_csv(args.output,sep="\t")
  print("Done.")
  ```

  + subset_by_col.py

  ```bash
  # single Normal
  python ../scripts/subset_by_col.py -i group_average.txt -o Normal.txt -w N -u NC
  # single Tumor
  python ../scripts/subset_by_col.py -i group_average.txt -o Tumor.txt -w T
  # single plasma, CRC-2399129 low mapping ratio after repair STAR mapping
  python ../scripts/subset_by_col.py -i group_average.txt -o plasma.txt -w plasma -u 2399129
  # all Average
  # Geneid|Length   Average_NC      Average_NC_PKU  Average_NC_ShH  Average_NC_ChQ  Average_CRC     Average_N       Average_T
  python ../scripts/subset_by_col.py -i group_average.txt -o Average.txt -w Average -u Unwanted
  
  ##make rnk
  cat Average.txt | awk '($2+0!=0){print $1,$6/$2}' > plasma_average_CRC_vs_average_HD.rnk
  cat Average.txt | awk '($3+0!=0){print $1,$6/$3}' > plasma_average_CRC_vs_average_HD_PKU.rnk
  cat Average.txt | awk '($7+0!=0){print $1,$8/$7}' > Tissue_average_Tumor_vs_average_Normal.rnk
  
  (base) [taoyuhuan@cnode Rankset_1027]$ paste Average.txt plasma.txt | head -1
  Geneid|Length	Average_NC	Average_NC_PKU	Average_NC_ShH	Average_NC_ChQ	Average_CRC	Average_N	Average_T	Geneid|Length	CRC-2416785-plasma	CRC-2418277-plasma	CRC-2418488-plasma	CRC-2418503-plasma	CRC-2418658-plasma
  (base) [taoyuhuan@cnode Rankset_1027]$ paste Average.txt plasma.txt | awk '($2+0!=0){print $1,$10/$2}' > CRC-2416785-plasma_vs_average_HD.rnk
  (base) [taoyuhuan@cnode Rankset_1027]$ paste Average.txt plasma.txt | awk '($2+0!=0){print $1,$11/$2}' > CRC-2418277-plasma_vs_average_HD.rnk
  (base) [taoyuhuan@cnode Rankset_1027]$ paste Average.txt plasma.txt | awk '($2+0!=0){print $1,$12/$2}' > CRC-2418488-plasma_vs_average_HD.rnk
  (base) [taoyuhuan@cnode Rankset_1027]$ paste Average.txt plasma.txt | awk '($2+0!=0){print $1,$13/$2}' > CRC-2418503-plasma_vs_average_HD.rnk
  (base) [taoyuhuan@cnode Rankset_1027]$ paste Average.txt plasma.txt | awk '($2+0!=0){print $1,$14/$2}' > CRC-2418658-plasma_vs_average_HD.rnk
  (base) [taoyuhuan@cnode Rankset_1027]$ paste Average.txt plasma.txt | awk '($3+0!=0){print $1,$10/$3}' > CRC-2416785-plasma_vs_average_HD_PKU.rnk
  (base) [taoyuhuan@cnode Rankset_1027]$ paste Average.txt plasma.txt | awk '($3+0!=0){print $1,$11/$3}' > CRC-2418277-plasma_vs_average_HD_PKU.rnk
  (base) [taoyuhuan@cnode Rankset_1027]$ paste Average.txt plasma.txt | awk '($3+0!=0){print $1,$12/$3}' > CRC-2418488-plasma_vs_average_HD_PKU.rnk
  (base) [taoyuhuan@cnode Rankset_1027]$ paste Average.txt plasma.txt | awk '($3+0!=0){print $1,$13/$3}' > CRC-2418503-plasma_vs_average_HD_PKU.rnk
  (base) [taoyuhuan@cnode Rankset_1027]$ paste Average.txt plasma.txt | awk '($3+0!=0){print $1,$14/$3}' > CRC-2418658-plasma_vs_average_HD_PKU.rnk
  
  
  ```

  

# RNA Splicing

### A. [MISO](https://miso.readthedocs.io/en/fastmiso/#ways-of-running-miso-and-associated-file-formats)

##### 000. installation and test

```bash
conda create -n miso python=2.7
conda activate miso
#install
pip install misopy
#check module, samtools is needed
module_availability
#test
python -m unittest discover misopy

```

```python
# prepera test files if error with no such files
cp -r /home/taoyuhuan/tools/miniconda3/envs/miso/misopy/* /home/taoyuhuan/tools/miniconda3/envs/miso/lib/python2.7/site-packages/misopy/
### 修改部分sam_to_bam.py(wd:/home/taoyuhuan/tools/miniconda3/envs/miso/lib/python2.7/site-packages/misopy/sam_to_bam.py)
##1.把绝对路径加到samtools前 2.修改sorted_filenames = ”%s.sorted.bam“ %(bam_filename.split(".bam")[0]) 3. cmd = "/home/taoyuhuan//tools/miniconda3/envs/miso/bin/samtools sort %s -o %s"加入-o，以适应新版本samtools

# Sort
    print "Sorting BAM file..."
    sorted_filename = "%s.sorted.bam" %(bam_filename.split(".bam")[0])
    cmd = "/home/taoyuhuan//tools/miniconda3/envs/miso/bin/samtools sort %s -o %s" %(bam_filename,
                                  sorted_filename)
    print "  - Executing: %s" %(cmd)
    os.system(cmd)

    # Index
    final_filename = "%s.bam" %(sorted_filename)
    print "Indexing BAM..."
    cmd = "/home/taoyuhuan//tools/miniconda3/envs/miso/bin/samtools index %s" %(final_filename)
    print "  - Executing: %s" %(cmd)
    os.system(cmd)
```

##### 00. Preparation

```bash
# py2.7 environment with samtools and bedtools
conda activate miso
```

##### 01. Insert Length of paired end sequencing

```bash
#!/bin/bash

input="/home/taoyuhuan/liquid/STAD/05.mapping_fq"
output="/home/taoyuhuan/liquid/STAD/14.insert_size_by_miso/exon_1000"
reference="/home/taoyuhuan/reference/forMISO/exon_1000"

echo $(date +%F%n%T)
exon_utils --get-const-exons /home/taoyuhuan/reference/forMISO/gencode.v34.annotation.gff3 --min-exon-size 1000 --output-dir ${reference}

for file in $(ls ${input}|grep sortedByCoord.out.bam)
	do
	echo ${file}
	echo $(date +%F%n%T)
	pe_utils --compute-insert-len ${input}/${file} ${reference}/gencode.v34.annotation.min_1000.const_exons.gff  --output-dir ${output}
done
```

##### 02. summary insert length

```bash
# !/bin/bash
cd /home/taoyuhuan/liquid/STAD/14.insert_size_by_miso/exon_1000
ls *.insert_len > name.txt
head -n 1 -q *.insert_len > insert_len.txt
paste name.txt insert_len.txt > insert_len_summary.txt
rm name.txt
rm insert_len.txt
```

##### 03. index gff, compute PSI and summary results

```bash
# !/bin/bash
input="/home/taoyuhuan/liquid/STAD/05.mapping_fq"
output="/home/taoyuhuan/liquid/STAD/15.miso_paired_end"
reference="/home/taoyuhuan/reference/forMISO/hg19/indexed_SE"
#index gff from MISO annotated gff
#index_gff --index /home/taoyuhuan/reference/forMISO/hg19/ /home/taoyuhuan/reference/forMISO/hg19/indexed_SE/

#readlines of pre-processed insert_len_summary.txt
cat /home/taoyuhuan/liquid/STAD/14.insert_size_by_miso/exon_1000/insert_len_summary.txt | while read line
do
#get filename,mean length and sdv from var(line,str)
file=`echo ${line}|awk -F '#' '{print $1}'` 
file=`echo ${file%.*}`
mean=`echo ${line}|awk -F '[#,=]' '{print $3}'`
sdv=`echo ${line}|awk -F '[#,=]' '{print $5}'`
# index bam
#samtools index ${input}/${file} ${input}/${file}.bai
#run miso paired end mode
echo 'file='${file},'mean_length='${mean},'sdv='${sdv}
echo $(date +%F%n%T)
mkdir ${output}/${file}
miso --run ${reference} ${input}/${file} --output-dir ${output}/${file} --read-len 150 --paired-end ${mean} ${sdv} --settings-filename /home/taoyuhuan/tools/miniconda3/misopy/settings/miso_settings.txt

summarize_miso --summarize-samples ${output}/${file} ${output}/${file}/summary_output/

done
```

+ Remember to set miso_settings.txt

```
[data]
filter_results = True
min_event_reads = 20
strand = fr-unstranded

[cluster]
cluster_command = qsub

[sampler]
burn_in = 500
lag = 10
num_iters = 5000
num_chains = 6
num_processors = 4
```

##### 04. comparison between two samples and filter diff splicing events

```bash
# compare
compare_miso --compare-samples NC_PKU-27Aligned.sortedByCoord.out.bam STAD_2163710Aligned.sortedByCoord.out.bam comparison
# filter
filter_events --filter NC_PKU-27Aligned.sortedByCoord.out.bam_vs_STAD_2163710Aligned.sortedByCoord.out.bam.miso_bf --num-inc 1 --num-exc 1 --num-sum-inc-exc 10 --delta-psi 0.20 --bayes-factor 10 --output-dir filtered/
```

##### 05. Sashimi plot

```bash
sashimi_plot --plot-event "chr9:34241183:34242106:+@chr9:34249777:34249959:+@chr9:34250656:34250757:+" /home/taoyuhuan/reference/forMISO/hg19/indexed /home/taoyuhuan/liquid/STAD/scripts/splicing/sashimi_plot_settings.txt --output-dir /home/taoyuhuan/liquid/STAD/15.1.splicing_plots/sashimi_plots

```

+ Remember to set sashimi_plot_settings.txt

```
[data]
# directory where BAM files are
bam_prefix = /home/taoyuhuan/liquid/STAD/05.mapping_fq/
# directory where MISO output is
miso_prefix = /home/taoyuhuan/liquid/STAD/15.miso_paired_end

bam_files = [
    "NC_PKU-27Aligned.sortedByCoord.out.bam",
    "NC_PKU-28Aligned.sortedByCoord.out.bam",
    "STAD_2163710Aligned.sortedByCoord.out.bam",
    "STAD_2244628Aligned.sortedByCoord.out.bam"]

miso_files = [
    "NC_PKU-27Aligned.sortedByCoord.out.bam",
    "NC_PKU-28Aligned.sortedByCoord.out.bam",
    "STAD_2163710Aligned.sortedByCoord.out.bam",
    "STAD_2244628Aligned.sortedByCoord.out.bam"]

[plotting]
# Dimensions of figure to be plotted (in inches)
fig_width = 7
fig_height = 5 
# Factor to scale down introns and exons by
intron_scale = 30
exon_scale = 4
# Whether to use a log scale or not when plotting
logged = False 
font_size = 6

bar_posteriors = False

# Max y-axis
ymax = 150

# Axis tick marks
nyticks = 3
nxticks = 4

# Whether to show axis labels
show_ylabel = True
show_xlabel = True

# Whether to plot posterior distributions inferred by MISO
show_posteriors = True 

# Whether to plot the number of reads in each junction
number_junctions = True

resolution = .5
posterior_bins = 40
gene_posterior_ratio = 5

# List of colors for read denisites of each sample
colors = [
    "#CC0011",
    "#CC0011",
    "#FF8800",
    "#FF8800"]

# Number of mapped reads in each sample
# (Used to normalize the read density for RPKM calculation)
coverages = [
    6830944,
    14039751,
    4449737, 
    6720151]

# Bar color for Bayes factor distribution
# plots (--plot-bf-dist)
# Paint them blue
bar_color = "b"

# Bayes factors thresholds to use for --plot-bf-dist
bf_thresholds = [0, 1, 2, 5, 10, 20]

##
## Names of colors for plotting
##
# "b" for blue
# "k" for black
# "r" for red
# "g" for green
#
# Hex colors are accepted too.
```

![image-20200702164237437](C:\Users\Tao\Documents\GitHub\data processing\sashimi_plot.png)

### B. [rMATS](http://rnaseq-mats.sourceforge.net/index.html)

```bash
sh run_rmats --b1 /home/taoyuhuan/liquid/STAD/scripts/splicing/rmats/NC2.txt --b2 /home/taoyuhuan/liquid/STAD/scripts/splicing/rmats/STAD2.txt --gtf /home/taoyuhuan/reference/forrMATS/gencode.v34.annotation.gtf -t paired --readLength 150 --nthread 4 --od /home/taoyuhuan/liquid/STAD/16.rMATS_Alternative_Splicing --tmp /home/taoyuhuan/liquid/STAD/16.rMATS_Alternative_Splicing/tmp
```

# Metagenomics

### matrix2[biom](https://blog.csdn.net/woodcorpse/article/details/84678543?utm_medium=distribute.pc_aggpage_search_result.none-task-blog-2~all~baidu_landing_v2~default-1-84678543.nonecase&utm_term=%E9%82%A3%E4%BA%9B%E8%BD%AF%E4%BB%B6%E5%8F%AF%E4%BB%A5%E5%AE%9E%E7%8E%B0biom%E6%A0%BC%E5%BC%8F%E7%9A%84%E8%BD%AC%E6%8D%A2&spm=1000.2123.3001.4430)

```bash
https://www.omicsclass.com/article/1314

biom convert -i table.biom -o table.from_biom_w_taxonomy.txt --to-tsv --header-key taxonomy


```



### A. picrust2

```bash
# install by conda *
conda create -n picrust2 -c bioconda -c conda-forge picrust2
# install by git

```

### B. [HUMAnN 2.0](https://github.com/biobakery/biobakery/wiki/humann2)

+ Installation

```bash
# add biobakery to first channel
conda config --add channels biobakery
# create env for hummann2
conda create -n humann2 python=2.7
# install by conda
conda install -c biobakery humann2
conda install -c bioconda metaphlan2=2.6.0
# [export metaphlan scripts to PATH](not required if use conda)(https://github.com/biobakery/biobakery/wiki/metaphlan2)
# export PATH=$MDIR:$MDIR/utils/:$PATH
# downlaod reference data
humann2_databases --download chocophlan full /home/taoyuhuan/reference/forHUMAnN2/
humann2_databases --download uniref uniref90_diamond /home/taoyuhuan/reference/forHUMAnN2/
# download utility_mapping for rename and regroup
 humann2_databases --download utility_mapping full /home/taoyuhuan/reference/forHUMAnN2/
# test
humann2 --input demo.fastq --output demo_fastq

# but error
# CRITICAL ERROR: Can not call software version for metaphlan2.py
# because conda install -c biobakery humann2 installed a higher version metaphlan2==2.7
# if conda install -c bioconda metaphlan2=2.6.0 installed humann2-2.8.1-py27_0
# we should re-install conda install -c biobakery humann2
# bioconda::humann2-2.8.1-py27_0 --> biobakery::humann2-0.11.2-py27_0
```

+ Usage

``` bash
# config database folder
humann2_config --update database_folders protein /home/taoyuhuan/reference/forHUMAnN2/uniref
humann2_config --update database_folders nucleotide /home/taoyuhuan/reference/forHUMAnN2/chocophlan
# run: humann2 --input {file} --output {dir}
humann2 --input /BioII/lulab_b/jinyunfan/projects/exSEEK/exSeek-dev/output/pico-final/unmapped/CRC-2124325/circRNA_1.fastq.gz --output CRC-2124325
# attach features to name
humann2_rename_table --input CRC-2124325_genefamilies.tsv --output CRC-2124325_uniref90_genefamilies-names.tsv --names uniref90
# regroup and rename
humann2_regroup_table -i CRC-2124325_genefamilies.tsv -g uniref90_ko -o CRC-2124325_ko_genefamilies.tsv
humann2_rename_table --input CRC-2124325_ko_genefamilies.tsv --names kegg-orthology --output CRC-2124325_ko_genefamilies-names.tsv
# cp to one dir and summary
cp /home/taoyuhuan/liquid/metagenomics/humann2_results/*/*_pathabundance.tsv /home/taoyuhuan/liquid/metagenomics/humann2_results_summary/
humann2_join_tables --input /home/taoyuhuan/liquid/metagenomics/humann2_results_summary/ --output Pathway_Abundance_summary.tsv
# delete repeat 
sed 's/_Abundance//g' Pathway_Abundance_summary.tsv > test.tsv
# normalize to cpm or relative abundance(relab)
humann2_renorm_table --input test.tsv --units relab --output Pathway_Abundance_normalized.tsv
# test by KW H-test
humann2_associate -i Cancertypes_associate_pathway_abundance_partial.tsv -m CANCERtypes -l CANCERtypes -t categorical -o /home/taoyuhuan/liquid/metagenomics/humann2_summary/cancertypes_pathway_KW-H-test.txt -f 0.2
```

+ plot by humann2

  ```bash
  humann2_barplot --input Cancertypes_associate_pathway_abundance_partial.tsv --focal-feature PWY-7219 --focal-metadatum CANCERtypes --last-metadatum CANCERtypes -x -s sum --output PWY-7219.png
  ```

+ plot by graphlan

  ```bash
  # install in humann2 env (python2.7)
  conda instal graphlan
  # make anno and tree for plot
  /home/taoyuhuan/tools/miniconda3/envs/humann2/bin/export2graphlan/export2graphlan.py --skip_rows 1,2 -i /home/taoyuhuan/liquid/metagenomics/humann2_summary/buglist_summary.tsv --tree merged_abundance.tree.txt --annotation merged_abundance.annot.txt --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 --annotations 5,6 --external_annotations 7 --min_clade_size 1
  # merge
  /home/taoyuhuan/tools/miniconda3/envs/humann2/bin/graphlan_annotate.py --annot merged_abundance.annot.txt merged_abundance.tree.txt merged_abundance.xml
  # plot
  /home/taoyuhuan/tools/miniconda3/envs/humann2/bin/graphlan.py --dpi 300 merged_abundance.xml merged_abundance.png --external_legends
  ```

  

### C. QIIME2

+ Installation

```bash
wget https://data.qiime2.org/distro/core/qiime2-2020.6-py36-linux-conda.yml
conda env create -n qiime2-2020.6 --file qiime2-2020.6-py36-linux-conda.yml
conda activate qiime2-2020.6
```

+ Usage

```bash

```

### D.mmvec

+ Installation

  ```bash
  conda activate qiime2-2020.6
  pip install git+https://github.com/biocore/mmvec.git
  qiime dev refresh-cache
  ```
  
+ [data](https://github.com/knightlab-analyses/multiomic-cooccurrences/tree/master/data) download

+ demo

+ 

# Shallow DNA seq for SNA detetcion

### WisecondorX

##### 00. Install

```bash
conda create -n WisecondorX python=2.7
conda activate WisecondorX
conda install -f -c conda-forge -c bioconda wisecondorx
## check by
WisecondorX
## any import error can be fixed by correct the import grammar acorrding to specific version

```



##### 01.test



# Test

1. Binomial test

   ```R
   #Binomial_pvalue
   output_pvalue=c()
   for(i in 1:nrow(LUAD_df)){
     p <- binom.test(x=as.numeric(LUAD_df[i,'Postive.in.Current.smoker']),n=as.numeric(LUAD_df[i,'Total.Current.smoker']),p=as.numeric(LUAD_df[i,'Proportion.in.Never']),alternative="greater")
     output_pvalue=rbind(output_pvalue,as.data.frame(p$p.value))
   }
   rownames(output_pvalue) <- LUAD_df$锘縯axo
   output_pvalue <- cbind(LUAD_df,output_pvalue)
   write.csv(output_pvalue,"p.csv")
   ```

2. 2 Proportion test in R？

   ```R
   #2 proportion_pvalue
   output_pvalue_2prop=c()
   for(i in 1:nrow(LUAD_df)){
     p <- prop.test(c(as.numeric(LUAD_df[i,'Postive.in.Never.smoker']),as.numeric(LUAD_df[i,'Postive.in.Current.smoker'])),c(as.numeric(LUAD_df[i,'Total.Never.smoker']),as.numeric(LUAD_df[i,'Total.Current.smoker'])))
     output_pvalue_2prop=rbind(output_pvalue_2prop,as.data.frame(p$p.value))
   }
   rownames(output_pvalue_2prop) <- LUAD_df$锘縯axo
   head(output_pvalue_2prop)
   output_pvalue_2prop <- cbind(output_pvalue,output_pvalue_2prop)
   write.csv(output_pvalue_2prop,"p_2prop.csv")
   ```

3. 2 Proportion test in python3.7 √

   ```python
   # !/bin/python
   # -*- coding: utf-8 -*-
   
   import numpy as np
   import pandas as pd
   from statsmodels.stats.proportion import proportions_ztest
   from sys import argv
   script, filename = argv
   i = 1
   df = pd.read_csv(filename,header=0)
   output = open("pvalue.txt",'w')
   output.write(filename+"p_2proportion\n")
   output.close()
   for i in range(len(df)):
   	count = np.array([df.loc[i,'Positive.in.low'],df.loc[i,'Positive.in.high']])
   	nobs = np.array([df.loc[i,'Total.CEA.low'],df.loc[i,'Total.high']])
   	stat, pval = proportions_ztest(count, nobs)
   	output = open("pvalue.txt",'a')
   	output.write(str(pval)+"\n")
   	output.close()
   ```
   
4. hyper geometry test

   ```R
   #phyper test(hypergeometry test)
   #phyper(q,M,N-M,k)
   #q=picked number of a specific genus(N=genus exist in one group)
   #M=exist number of the genus in both group(M=genus exist in group A+B)
   #N=all samples
   #k=trails number(k>q,k=sample number in one group)
   ?phyper
   output_pvalue_phyper=c()
   for(i in 1:nrow(LUAD_df)) {
     if(as.numeric(LUAD_df[i,'Proportion.in.Current'])>as.numeric(LUAD_df[i,'Proportion.in.Never'])) {
       M=as.numeric(LUAD_df[i,'Postive.in.Current.smoker'])+as.numeric(LUAD_df[i,'Postive.in.Never.smoker'])
       p <- phyper(as.numeric(LUAD_df[i,'Postive.in.Current.smoker'])-1,M,34-M,as.numeric(LUAD_df[i,'Total.Current.smoker']),lower.tail=FALSE)
       output_pvalue_phyper=rbind(output_pvalue_phyper,as.data.frame(p))
       }
     else{
       M=as.numeric(LUAD_df[i,'Postive.in.Current.smoker'])+as.numeric(LUAD_df[i,'Postive.in.Never.smoker'])
       p <- phyper(as.numeric(LUAD_df[i,'Postive.in.Never.smoker'])-1,M,34-M,as.numeric(LUAD_df[i,'Total.Never.smoker']),lower.tail=FALSE)
       output_pvalue_phyper=rbind(output_pvalue_phyper,as.data.frame(p))
       }
   }
   rownames(output_pvalue_phyper) <- LUAD_df$锘縯axo
   head(output_pvalue_phyper)
   output_pvalue_phyper <- cbind(output_pvalue_2prop,output_pvalue_phyper)
   write.csv(output_pvalue_phyper,"p_phyper.csv")
   ```

5. FDR_BH

   ```R
   #BH_fdr
   pvalue_to_fdr <- read.csv("LUAD_pvalue_filtered.csv",header=T)
   pvalue_to_fdr <- as.numeric(unlist(pvalue_to_fdr$p.p.value))
   length(pvalue_to_fdr)
   
   fdr <- p.adjust(pvalue_to_fdr, method = "BH", n = length(pvalue_to_fdr))
   as.data.frame(fdr)
   pvalue_to_fdr <- cbind(pvalue_to_fdr,fdr)
   write.csv(pvalue_to_fdr,"LUAD_filtered.csv")
   
   length(pvalue_to_fdr$p.p.value)
   ```



# Tips

### 00. Update R

```R
# 安装包"installr"
install.packages("installr")
# 导入包
library(installr)
# 升级
updateR()
# 包升级
update.packages()
# 选择镜像
options(repos=structure(c(CRAN="https://cran.cnr.berkeley.edu/")))
```

### 01. Conda

```bash
##### channels
# add channel
conda config --add channels biobakery
#change channel
pip install xxx -i https://pypi.tuna.tsinghua.edu.cn/simple
清华：https://pypi.tuna.tsinghua.edu.cn/simple
阿里：http://mirrors.aliyun.com/pypi/simple/
中国科技大学 https://pypi.mirrors.ustc.edu.cn/simple/
华中理工大学：http://pypi.hustunique.com/
山东理工大学：http://pypi.sdutlinux.org/ 
豆瓣：http://pypi.douban.com/simple/

###### envs
# check envs in conda
conda envs list
conda info --envs
# create python2.7 environment
conda create --name py2 python==2.7
# enter envs
conda activate py2
# exit envs
conda deactivate
# solve no name for environment（conflicted path with conda）
conda create -n rmats --clone /home/taoyuhuan/tools/rMATS/rmats-turbo/conda_envs/rmats
# remove environment with only path
conda remove -p /home/taoyuhuan/tools/rMATS/rmats-turbo/conda_envs/rmats --all
# check packages in py2 envs
pip list
# install packages in py2 envs
pip install xxx
# freeze envs packages and verison, preserve in a .txt, and re-install
pip freeze > py2.txt
pip install -r py2.txt
```

### 02. environment path

```bash
#check environment path
echo $PATH
# check the path of software in-use
which metaphlan.py
```

### 03. check system

```bash
# check system
lscpu
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                16
On-line CPU(s) list:   0-15
Thread(s) per core:    2
Core(s) per socket:    4
Socket(s):             2
NUMA node(s):          2
Vendor ID:             GenuineIntel
CPU family:            6
Model:                 79
Model name:            Intel(R) Xeon(R) CPU E5-2637 v4 @ 3.50GHz
Stepping:              1
CPU MHz:               1201.074
CPU max MHz:           3700.0000
CPU min MHz:           1200.0000
BogoMIPS:              6999.38
Virtualization:        VT-x
L1d cache:             32K
L1i cache:             32K
L2 cache:              256K
L3 cache:              15360K
NUMA node0 CPU(s):     0,2,4,6,8,10,12,14
NUMA node1 CPU(s):     1,3,5,7,9,11,13,15
Flags:                 fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc aperfmperf eagerfpu pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 fma cx16 xtpr pdcm pcid dca sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand lahf_lm abm 3dnowprefetch epb cat_l3 cdp_l3 invpcid_single intel_pt spec_ctrl ibpb_support tpr_shadow vnmi flexpriority ept vpid fsgsbase tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm cqm rdt_a rdseed adx smap xsaveopt cqm_llc cqm_occup_llc cqm_mbm_total cqm_mbm_local dtherm ida arat pln pts
# check real time cpu (live update) 
top

# check file/dir size
du -h --max-depth=1 /BioII/lulab_b/jinyunfan/projects/exSEEK/exSeek-dev/output/pico-final/unmapped/*/circRNA*
```

### 04. reference download

+ ensemble [GFF]([ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/](ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/))
+ encode [GFF](https://www.gencodegenes.org/human/)
+ UCSC GTF FASTA
+ NCBI 

### 05. log

```bash
nohup command > name.log 2>&1 &
```

### 06. R loop

```R
setwd("C:/Users/Tao/Desktop/")
LUAD_df <- read.csv("LUAD_df.csv",header=TRUE)
head(LUAD_df)

output_pvalue=c()
for(i in 1:nrow(LUAD_df)){
  p <- binom.test(x=as.numeric(LUAD_df[i,'Postive.in.Current.smoker']),n=as.numeric(LUAD_df[i,'Total.Current.smoker']),p=as.numeric(LUAD_df[i,'Proportion.in.Never']),alternative="greater")
  output_pvalue=rbind(output_pvalue,as.data.frame(p$p.value))
}

rownames(output_pvalue) <- LUAD_df$锘縯axo
write.csv(output_pvalue,"p.csv")
```

### 07. Bus error

```bash
##Problem: import numpy then exit python, report Bus error; conda env list, report bus error;
##solution:
#check storage
df -h
#or check tmp
df -h /tmp/
# if 100% used, delete something
```

### 08. Count genetype in GTF

```bash
cat gencode.gtf |perl -alne '{next unless $F[2] eq "gene" ;/gene_type "(.*?)";/;print $1}'  |sort |uniq -c > genetype.txt

##output
gencode.gtf
    14 IG_C_gene
      9 IG_C_pseudogene
     37 IG_D_gene
     18 IG_J_gene
      3 IG_J_pseudogene
      1 IG_pseudogene
    144 IG_V_gene
    188 IG_V_pseudogene
  16899 lncRNA
   1881 miRNA
   2212 misc_RNA
      2 Mt_rRNA
     22 Mt_tRNA
     49 polymorphic_pseudogene
  10169 processed_pseudogene
  19954 protein_coding
     18 pseudogene
      8 ribozyme
     47 rRNA
    497 rRNA_pseudogene
     49 scaRNA
      1 scRNA
    943 snoRNA
   1901 snRNA
      5 sRNA
   1058 TEC
    500 transcribed_processed_pseudogene
    138 transcribed_unitary_pseudogene
    941 transcribed_unprocessed_pseudogene
      2 translated_processed_pseudogene
      1 translated_unprocessed_pseudogene
      6 TR_C_gene
      4 TR_D_gene
     79 TR_J_gene
      4 TR_J_pseudogene
    106 TR_V_gene
     33 TR_V_pseudogene
     97 unitary_pseudogene
   2615 unprocessed_pseudogene
      1 vault_RNA

```

