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

