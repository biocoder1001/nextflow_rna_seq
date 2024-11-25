# nextflow_rna_seq

This repository contains a Nextflow pipeline designed for RNA-Seq analysis. The pipeline is tailored to run on a SLURM cluster and automates the complete workflow from raw data preprocessing to differential expression and downstream PCA analysis.

Overview
The pipeline processes FASTQ files, performs essential preprocessing steps, and conducts differential expression analysis for a subset of samples specified by the user. It is capable of handling multi-level factor analysis to explore differential gene expression comprehensively.

Features
Preprocessing: Automates quality control, alignment, and other preparatory steps for all input FASTQ files.
Flexible Sample Selection: Only the samples specified in the samples.txt file are included in the differential expression and PCA analyses.
Customizable Parameters: Parameters can be set in the helper file or overridden via the command line.
Multi-Factor Analysis: Supports multi-level factor analysis for robust differential expression studies.

Input Files
Required Inputs:
FASTQ Files: Raw RNA-Seq data files to be processed.
samples.txt: A tab-delimited file specifying:
Sample ID: Unique identifier for each sample.
Donor: Indicates the donor source for grouping or batch correction.
Condition: Experimental condition (e.g., treatment, control, etc.).
Example of samples.txt:
Sample ID	Donor	Condition
Sample_001	Donor1	Control
Sample_002	Donor1	Treated
Sample_003	Donor2	Control

Workflow
Preprocessing:

Processes all provided FASTQ files.
Includes steps  trimming, alignment, and quantification.

Subset Selection:
Filters the samples for differential expression based on the samples.txt file.

Differential Expression Analysis:
Identifies differentially expressed genes between specified conditions.
Includes support for multi-factor analysis.

PCA Analysis:
Performs Principal Component Analysis on the filtered dataset.

Customization
Command-Line Overrides: All parameters can be overridden directly through the command line, allowing for dynamic configuration.
 

