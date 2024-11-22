#!/usr/bin/env nextflow
/*
* This is a rnaseq pipeline
*/

nextflow.enable.dsl=2

params.index = '/fsimb/groups-external/sfbcardosogr/ishita/ensembl/grch38/canonical/index/star/2.7.3a'
params.reads = '/fsimb/groups-external/sfbcardosogr/ishita/maruthi_data/fastq_files/*.fastq.gz'
params.gtf_file='/fsimb/groups-external/sfbcardosogr/ishita/ensembl/grch38/canonical/annotation/Homo_sapiens.GRCh38.98.gtf'
params.samples='/fsimb/groups-external/sfbcardosogr/ishita/maruthi_data/samples.txt'
params.basedir = '/fsimb/groups-external/sfbcardosogr/ishita/maruthi_data'
log.info """\
	RNASEQ-NF-PIPELINE
	==================
	ISHITA | CARDOSO LAB | TECHNICAL UNIVERSITY OF DARMSTADT
	"""

// Define channels
reads_ch = Channel.fromPath(params.reads,checkIfExists: true)
index_ch = Channel.value(params.index)
gtf_ch = Channel.value(params.gtf_file)
sample_selection = Channel.value(params.samples)

// Process definition
process fastqc_qc_alignment {
    publishDir  "${params.basedir}/fastc_results", mode: 'copy'
    input:
    path reads

    output:
    path "fastqc_files" 

    script:
    """
    # Create directories for output
    echo '${reads}'
    mkdir -p fastqc_files

    # Run FastQC on reads
    fastqc ${reads} -o fastqc_files
    """
}

process cutadapt {
	publishDir  "${params.basedir}/cutadapt_results", mode: 'copy'
	input:
	path reads
	
	output:
	path "cutadapt/*_cutadapt.fastq.gz"

	script:
	"""
	mkdir -p cutadapt
	cutadapt --adapter Illumina=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --overlap=5 --minimum-length 30 --nextseq-trim=20 --error-rate 0.1  --output=cutadapt/${reads.baseName}.cutadapt.fastq.gz $reads 
	"""
}

process alignment_quant {
	publishDir  "${params.basedir}/aligned_bam_files", mode: 'copy'
	input:
	path reads
	path index
	path gtf
		
	output:
	path "bam_files/*.bam", emit: bam_files
	path "count_files/*.ReadsPerGene.out.tab", emit: count_matrix 
	
	script:
	"""
	mkdir -p bam_files
	mkdir -p count_files
	STAR --runThreadN 4 \
	 --quantMode GeneCounts\
	 --runMode alignReads\
	 --outFilterMismatchNmax 999\
	 --sjdbOverhang 49\
	 --sjdbGTFfile $gtf\
         --genomeDir $index \
         --readFilesIn $reads \
         --readFilesCommand zcat \
         --outFileNamePrefix bam_files/${reads.baseName}. \
         --outSAMtype BAM SortedByCoordinate
	 mv bam_files/*.ReadsPerGene.out.tab count_files/
	"""
}

process counts_mat {
	publishDir  "${params.basedir}/coutn_files", mode: 'copy'
	input:
	path counts

	output:
	path "count_files/*.tsv"

	script:
	"""
	mkdir -p count_files
	sed -i '1,4d' ${counts}
	cut -f1,2 ${counts} > ${counts.baseName}.tsv
	mv ${counts.baseName}.tsv count_files/ 
	"""
}

process sample_file_extract {
	publishDir  "${params.basedir}/filtered_counts", mode: 'copy'
	input:
	path count_files
	path samples
	
	output:
	path  "*.tsv"

	script:
    	"""
	sample_list=\$(cut -f1 ${samples})
	sample_base=\$(basename ${count_files} | cut -d'.' -f1)
    	found=false
    	while IFS= read -r sample|| [[ -n "\$sample" ]]; do
        	sample=\$(echo "\$sample" | tr -d '\r')
        	if [[ "\$sample_base" == "\$sample" ]]; then
            		echo "Match found for \$sample_base"
            # No need to copy since the file is already in the working directory
            		found=true
			cp ${count_files} \${sample_base}.tsv
            		break
        	fi
    	done <<< \${sample_list}

    	if [[ "\$found" == "false" ]]; then
       		 echo "No match found for \$sample_base in samples file"
		 touch empty.tsv
    	fi	
    	"""
}


process DeSeq {
	publishDir  "${params.basedir}/combined_counts", mode: 'copy'
	input:
	val count_file_list 
	path samples
	
	output:
	path "combined_count_matrix.tsv"

	script:
	"""
	#!/usr/bin/env Rscript
	library(tools)
	sampleFiles <- c(${count_file_list.join('","').replaceAll('^', '"').replaceAll('$', '"')})
	cond <- read.delim("$samples",  header = TRUE)
	sampleTable <- data.frame(sample = cond[,1], file = sampleFiles, group = cond[,2], donor = cond[,3])
	write.table(sampleTable, file="combined_count_matrix.tsv", sep="\t", row.names=FALSE, quote=FALSE)
	"""
}
	
process deseq_i {
	publishDir  "${params.basedir}/counts_file_assay", mode: 'copy'
	input:
	path metadata
	
	output:
	'*.tsv'

	script:
	"""
	#!/usr/bin/env Rscript
	library(DESeq2)
	targets <- read.delim("$metadata",  header = TRUE)
	print(targets)
	directory = '/fsimb/groups-external/sfbcardosogr/ishita/maruthi_data/filtered_counts'
	group <- factor(targets\$group)
	donor <- factor(targets\$donor)
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = targets,
				       directory = directory,
				       design = ~group)
	dds <- DESeq(ddsHTSeq)
	res <- results(dds, alpha=0.01)
	write.table(res, "deg.tsv", sep="\t", row.names=FALSE, quote=FALSE) 
	"""
}

// Workflow definition
workflow {
    // Execute the process and capture the output into a channel
    fastqc_qc_alignment(reads_ch)
    trimmed_ch=cutadapt(reads_ch)
    alignment_ch=alignment_quant(trimmed_ch,index_ch,gtf_ch)
    processed_counts = counts_mat(alignment_ch.count_matrix)
    result=sample_file_extract(processed_counts,sample_selection).filter{file ->file.name =~ /D.*\.tsv/}
    count_files_list = result.collect().name
    count_files_list.view()
    deseq_ch=DeSeq(count_files_list,sample_selection)
    deseq_analysis = deseq_i(deseq_ch)
    
}

workflow.onComplete {
	println (workflow.success ? "Done! Wrapping up..." : "Oooops!!! Somehing just went wrong! Go check it")
}
