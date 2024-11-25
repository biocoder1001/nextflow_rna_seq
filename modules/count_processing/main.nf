process counts_mat_processing {
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

