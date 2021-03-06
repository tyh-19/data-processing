# **RNA different expression**

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

# RNA Splicing

### A. [MISO](https://miso.readthedocs.io/en/fastmiso/#ways-of-running-miso-and-associated-file-formats)

##### 00. Preparation

```bash
# py2.7 environment with samtools and bedtools
conda activate py2
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
```



### 03. reference download

+ ensemble [GFF]([ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/](ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/))
+ encode [GFF](https://www.gencodegenes.org/human/)
+ UCSC GTF FASTA
+ NCBI 

