process fastqc {
    publishDir  "${params.basedir}/fastc_results", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_files"

    script:
    def single = reads instanceof Path
    if (!single)	
    """
    mkdir -p fastqc_files
    fastqc ${reads[0]} ${reads[1]} -o fastqc_files
    """
    else 
    """
    fastqc ${reads} -o fastqc_files
    """
}

