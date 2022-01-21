#! /usr/bin/env Rscript

input_dir <- "/data/taoyuhuan/projects/SNP_DNA_RNA/extract_VEP/VEP_RNA"
output_dir <- "/data/taoyuhuan/projects/SNP_DNA_RNA/extract_VEP/concordant_VEP"
input_concordance <- "/data/taoyuhuan/projects/SNP_DNA_RNA/individual_compare"
samples <- as.character(lapply(strsplit(dir(input_dir),".",fixed = TRUE), function(x) x[1]))
samples <- gsub("-pico","",samples)
samples <- gsub("-wgs","",samples)
j=1
while(j<=length(samples)){
  print(samples[j])
  VEP <- read.table(paste0(input_dir,"/",samples[j],"-pico.rmEDIT.filtered.VEP.txt"),sep = "\t")
  concordent <- read.table(paste0(input_concordance,"/",samples[j],".concordant.vcf.diff.sites_in_files"),sep = "\t",header = TRUE)
  concordent <- concordent[which(concordent$IN_FILE=="B"),]
  #SNPs <- paste0(concordent$CHROM,"_",concordent$POS1,"_",concordent$REF1,"/",concordent$ALT1)
  SNPs <- paste0(concordent$CHROM,":",concordent$POS1)
    
  i=1
  concordent_VEP={}
  while(i<=length(SNPs)){
  SNPs_VEP <- VEP[grep(SNPs[i],VEP$V2),]
  concordent_VEP <- rbind(concordent_VEP,SNPs_VEP)
  i=i+1
  } 
  colnames(concordent_VEP) <- c("Uploaded_variation","Location","Allele","Gene","Feature","Feature_type","Consequence","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","IMPACT","DISTANCE","STRAND","FLAGS","VARIANT_CLASS","SYMBOL","SYMBOL_SOURCE","HGNC_ID","BIOTYPE","CANONICAL","MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM","GENE_PHENO","SIFT","PolyPhen","EXON","INTRON","DOMAINS","miRNA","HGVSc","HGVSp","HGVS_OFFSET","AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF","gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF","MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS")
  write.table(concordent_VEP,paste0(output_dir,"/",samples[j],".concordant.VEP.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
j=j+1
}
