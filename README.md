# Data-processing
- scripts of RNAseq data

# Multiomics
- **scripts of multiomics RNA paper (Figure 1 - Figure 5)**
  - Figure 1-5
    - Figure1.Sampel_Summary.R
    - Figure2.Omics_correlation_coverage.R
    - Figure3.Multiple alteration profiling.R
    - Figure4.Differential analysis of multiple alterations.R
    - Figure5.Alteration combination.R```

  - Further exploration
    - multiomics.R

# RP-mRNA
- scripts of RP-mRNA related analysis
  - exploration scripts
    - selectFlows.R

  - Orignazied scripts
    - RP-mRNA.R

# Useful scripts
- 00.Genome_annotation
- 01.RNA_features
  - APA: alternative polyadenylation
  - ASE: allelic specific expression
  - Alternative promoter
  - chimeric RNA
  - Editing
  - Insert length
  - Microbe
  - RNA SNP
  - RNA splicing
  - RNA transposon elements
  - Normailzed Bigwig
- 02.Subset_Norm
- 03.Differential_Expression
- 04.GeneCoverage
- 05.Correlation
- 06.MachineLearning
  - LOO_v3_yuhuan.R: Random forest leave-one-out cross validation.
  - Combination.R: combine different alteration types in model.
  - Sensitivity_at0.95speci.R: calculate sensitivity at 97% specificity
  - `*.sh: scripts on cluster 
- 07.Origin
  - Origin.R: exploration scripts
  - Signature_gene_matrix_v2.R: generate sig.gene ref from given expression matrix and plot
  - Multiomics_paired.R, exoRBase.R, GSE*.R: deconvolution scripts
  - Cibersort.R: source code for deconvolution
- Plots
  - Plot scripts, such as Dot heatmap
- 2021_simple scripts
  - 2021 simple scripts
- 2020_old pipelines
  - 2020 old pipelines
- CI
  - Confidence Interval calculate tools
- PCAWGplot
  - Data prepare for PCAWG plot: detection ratio for 95% outliers(relative to healthy group) in cancer group.
- count_matrix_summarization
  - Summarize TCGA counts
- qPCR tools
  - A uniqform plot tool for qPCR Ct value comparison

# small_cfRNAseq
- pipeline for small cfRNAseq
