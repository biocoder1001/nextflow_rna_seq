process fastqc {
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

