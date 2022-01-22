# Data-processing
- **Scripts of RNAseq data**

# Multiomics
- **Scripts of multiomics RNA paper (Figure 1 - Figure 5)**
  - Figure 1-5
    - Figure1.Sampel_Summary.R
    - Figure2.Omics_correlation_coverage.R
    - Figure3.Multiple alteration profiling.R
    - Figure4.Differential analysis of multiple alterations.R
    - Figure5.Alteration combination.R

  - Further exploration
    - multiomics.R

# RP-mRNA
- **Scripts of RP-mRNA related analysis**
  - Exploration scripts
    - selectFlows.R

  - Orignazied scripts
    - RP-mRNA.R

# Useful scripts
- **00.Genome_annotation**
- **1.RNA_features**
  - APA: alternative polyadenylation
  - ASE: allelic specific expression
  - Alternative promoter
  - Chimeric RNA
  - Editing
  - Insert length
  - Microbe
  - RNA SNP
  - RNA splicing
  - RNA transposon elements
  - Normailzed Bigwig
- **02.Subset_Norm**
- **03.Differential_Expression**
- **04.GeneCoverage**
- **05.Correlation**
- **06.MachineLearning**
  - LOO_v3_yuhuan.R: Random forest leave-one-out cross validation.
  ```
  usage: LOO_v3_yuhuan.R [-h] -i MATRIX -f FEATURENUMBER -p {CRC,STAD,NC} -n
                       {CRC,STAD,NC} -o OUTDIR
  Random Forest LOO
  optional arguments:
  -h, --help            show this help message and exit
  -i MATRIX, --matrix MATRIX
                        input count matrix. Rows are genes. Columns are
                        samples.
  -f FEATURENUMBER, --featureNumber FEATURENUMBER
                        feature number in final model. 200 is suggested. If
                        total features are less than 200, use all features.
  -p {CRC,STAD,NC}, --positive {CRC,STAD,NC}
  -n {CRC,STAD,NC}, --negative {CRC,STAD,NC}
  -o OUTDIR, --outdir OUTDIR
  ```
  - Combination.R: combine different alteration types in model.
  ```
  usage: Combination.R [-h] -i INPUT_DIR -f ALTERATIONNUMBER -p {CRC,STAD,NC} -n
                     {CRC,STAD,NC} -o OUTPUT_DIR
  Merge model
  optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
                        input directory, which contains results like
                        ./Expression/STADvsNC/ from LOO_v3_yuhuan.R.
  -f ALTERATIONNUMBER, --AlterationNumber ALTERATIONNUMBER
                        Alteration number in final model.
  -p {CRC,STAD,NC}, --positive {CRC,STAD,NC}
  -n {CRC,STAD,NC}, --negative {CRC,STAD,NC}
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
  ```
  - Sensitivity_at0.95speci.R: calculate sensitivity at 97% specificity
  ```
  usage: Sensitivity_at0.95speci.R [-h] -i INPUT_DIR -p {CRC,STAD,NC} -n
                                 {CRC,STAD,NC} -o OUTPUT_DIR

  optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
                        input dir contains a dataframe named
                        prob_final_LOO.txt. colnames are positive and negative
                        group. rownames are samples which named
                        "positive.xxx"/"negative.xxx".
  -p {CRC,STAD,NC}, --positive {CRC,STAD,NC}
  -n {CRC,STAD,NC}, --negative {CRC,STAD,NC}
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
  ```
  - *.sh: scripts on cluster 
- **07.Origin**
  - Origin.R: exploration scripts
  - Signature_gene_matrix_v2.R: generate sig.gene ref from given expression matrix and plot
  ```
  usage: Signature_gene_matrix_v2.R [-h] -m MATRIX -s SAMPLE_INFO
                                  [-e TPM_CUTOFF] [-ts TSS_CUTOFF] [-n NUMBER]
                                  -o OUTDIR
  Signature Matrix v2
  optional arguments:
  -h, --help            show this help message and exit
  -m MATRIX, --matrix MATRIX
                        Input count matrix. Rows are genes(TPM). Columns are
                        samples.
  -s SAMPLE_INFO, --sample_info SAMPLE_INFO
                        Sample information for each sample. CSV file with 2
                        columns. Column 1: sample, Column 2: group
  -e TPM_CUTOFF, --TPM_cutoff TPM_CUTOFF
                        TPM cutoff, maximum expressionn larger than this
                        cutoff. Default: 1.
  -ts TSS_CUTOFF, --TSS_cutoff TSS_CUTOFF
                        TSS cutoff, minimum tissue specific score larger than
                        this cutoff. Default: 0.5.
  -n NUMBER, --Number NUMBER
                        Signature gene number. Default: 20.
  -o OUTDIR, --outdir OUTDIR
                        output directory
  ```
  - Multiomics_paired.R, exoRBase.R, GSE*.R: deconvolution scripts
  - Cibersort.R: source code for deconvolution
- **CI**
  - Confidence Interval calculate tools
  - Confidence_Interval.R
  ```
  Usage: Confidence_Interval.R [options]
  Options:
	-i INPUT, --input=INPUT
		Input file names 'miso_Insert_length_summary.txt' or other file contain a string of numbers each column, split by tab
	-o OUTPUT, --output=OUTPUT
		Output directory and file name
	-b BOOTSTRAP_NUM, --bootstrap_num=BOOTSTRAP_NUM
		Bootstrap number, default is 1000
	-c CONFIDENCE_LEVEL, --confidence_level=CONFIDENCE_LEVEL
		Confidence level, defaut is 0.95
	-h, --help
		Show this help message and exit
  ```
- **PCAWGplot**
  - Data prepare for PCAWG plot: detection ratio for 95% outliers(relative to healthy group) in cancer group.
  - FDR.R
  ```
  usage: FDR.R [-h] -i MATRIX -o OUTDIR
  Outlier detection
  optional arguments:
  -h, --help            show this help message and exit
  -i MATRIX, --matrix MATRIX
                        input count matrix. Rows are features. Columns are
                        samples. Tab seperated.
  -o OUTDIR, --outdir OUTDIR
                        output directory(full path end with /). Make sure
                        directory exists.
  ```
- **count_matrix_summarization**
  - Summarize TCGA counts
- **Plots**
  - Plot scripts, such as Dot heatmap
- **qPCR tools**
  - A uniqform plot tool for qPCR Ct value comparison
- **2021_simple scripts**
  - 2021 simple scripts
- **2020_old pipelines**
  - 2020 old pipelines

# small_cfRNAseq
- **pipeline for small cfRNAseq**
