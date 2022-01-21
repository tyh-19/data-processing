```R
##waterfall plot by maftools
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("maftools")
  
  library(maftools)
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/SNP/annotation_filteredPBMC/")
  var.annovar.maf = annovarToMaf(annovar = "CRC_sample_CRC_COSMIC_cancergenes_filtered.txt", 
                                 Center = 'NA', 
                                 refBuild = 'hg38', 
                                 tsbCol = 'Tumor_Sample_Barcode', 
                                 table = 'refGene',
                                 sep = "\t")
  write.table(var.annovar.maf,file="var_annovar.maf",quote= F,sep="\t",row.names=F)
  var_maf = read.maf(maf ="var_annovar.maf")
  
  var_maf@clinical.data$Type <- as.character(lapply(strsplit(as.character(var_maf@clinical.data$Tumor_Sample_Barcode),"-"),function(x) x[1]))
  
  var_maf@clinical.data
  
  plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median')
  
  oncoplot(maf = var_maf, top = 25 ,showTumorSampleBarcodes = F,clinicalFeatures = "Type",colbar_pathway = FALSE,
           sortByAnnotation = TRUE,draw_titv = TRUE)
  
  laml.titv = titv(maf = var_maf, plot = FALSE, useSyn = TRUE)
  plotTiTv(res = laml.titv)
  
  somaticInteractions(maf = var_maf, top = 10, pvalue = c(0.05, 0.1))
  
  #lollipop plot for APC
  lollipopPlot(
    maf = var_maf,
    gene = 'APC',
    AACol = 'aaChange',
    showMutationRate = TRUE,
    #labelPos = "all",
    #refSeqID = "NM_000546",
    printCount = TRUE,
    showDomainLabel = TRUE
  )
```

