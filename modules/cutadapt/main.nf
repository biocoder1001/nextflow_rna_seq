process cutadapt {
        publishDir  "${params.basedir}/cutadapt_results", mode: 'copy'
        input:
        tuple val(sample_id), path(reads)

        output:
        path "cutadapt/*_cutadapt.fastq.gz"

        script:
        def single = reads instanceof Path
	if (!single)
        """
        mkdir -p cutadapt
        cutadapt --adapter Illumina=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --overlap=5 --minimum-length 30 --nextseq-trim=20 --error-rate 0.1 -o cutadapt/${reads[0].baseName}_cutadapt.fastq.gz  -p cutadapt/${reads[1].baseName}_cutadapt.fastq.gz ${reads[0]} ${reads[1]}
        """

	else
	"""
	mkdir -p cutadapt
	cutadapt --adapter Illumina=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --overlap=5 --minimum-length 30 --nextseq-trim=20 --error-rate 0.1 --output=cutadapt/${reads.baseName}_cutadapt.fastq.gz  ${reads} 
	"""
}

