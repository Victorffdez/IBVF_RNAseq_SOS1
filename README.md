# RNASeq pipeline. From FASTQ to DEGs and GO analysis.
This repository contains all the code used in the analysis of RNA-seq data from rice (_O. sativa_). You can see two different directories:

- **Pre-processing, Alignment & Quantification.** This folder contains the Bash scripts corresponding to the different stages of RNA sequencing data processing using the new Tuxedo protocol: pre-processing (*RawData_Fastqc.sh* & *TrimmedData_Fastqc.sh*), alignment (*Alignment_hisat2.sh*), assembly and quantification (*Assembly_&_Quantification_stringtie.sh*).
  
- **Differential Expression Analysis & GO Enrichment.**  R scripts to perform all the steps necessary to analyze the gene expression matrix obtained in the previous steps. Here, among many options, you will find scripts to normalize the raw data, perform an exhaustive exploratory analysis, use _limma_ to obtain DEGs between conditions and perform a GO term enrichment analysis.


Please read all the comments in the scripts carefully. If you have any questions, please do not hesitate to contact me at the follow e-mail address:

**Víctor J. Fernández Ramírez** 

*victor.fernandez@imibic.org* 

*Bioinformatician* 
