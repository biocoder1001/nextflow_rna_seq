process metadata_prep {
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


