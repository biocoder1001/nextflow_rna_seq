# nextflow_rna_seq

This is a nextflow pipeline for the rna seq analysis, it currently runs on the slurm cluster. It takes the fastq files, does all the 
required preprocessing, followed by only taking the samples that are required for the differential expression analysis and further 
downstream PCA analysis. The helper file explains the parameters required for the file. Although these params can be overwritten if you give it in the 
command line. 

