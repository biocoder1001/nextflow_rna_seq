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
	if ("donor" %in% colnames(cond)) {
    	    sampleTable <- data.frame(sample = cond[,1], file = sampleFiles, group = cond[,3], donor = cond[,2])
	} else {
	     sampleTable <- data.frame(sample = cond[,1], file = sampleFiles, group = cond[,2])
	   }
        write.table(sampleTable, file="combined_count_matrix.tsv", sep="\t", row.names=FALSE, quote=FALSE)
        """
}


