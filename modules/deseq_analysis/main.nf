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
	Rscript ${params.basedir}/DeSeq.R ${metadata} ${params.basedir} ${ref}  ${group_1} ${group_2}  downregulated_genes.tsv upregulated_genes.tsv
        """
}

