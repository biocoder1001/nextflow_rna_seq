process deseq_analysis {
	publishDir  "${params.basedir}/counts_file_assay", mode: 'copy'
        input:
        path metadata
	val ref
	val group_1
	val group_2
        
	output:
        path "downregulated_genes.tsv"
	path "upregulated_genes.tsv"
	path "DE_DESeq2.pdf"
	
        script:
        """
	Rscript ${params.basedir}/deseq_edited.R -m ${metadata} -d ${params.basedir} -r ${ref}  -f ${group_1} -s ${group_2}  -l downregulated_genes.tsv -u upregulated_genes.tsv
        """
}

