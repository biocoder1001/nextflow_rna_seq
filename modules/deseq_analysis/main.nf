process deseq_analysis {
        publishDir  "${params.basedir}/counts_file_assay", mode: 'copy'
        input:
        path metadata

        output:
        '*.tsv'

        script:
        """
	Rscript ${params.basedir}/DeSeq.R ${metadata}
        """
}

