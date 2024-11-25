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

