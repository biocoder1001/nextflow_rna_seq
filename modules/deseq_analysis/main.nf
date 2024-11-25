process deseq_analysis {
        publishDir  "${params.basedir}/counts_file_assay", mode: 'copy'
        input:
        path metadata

        output:
        '*.tsv'

        script:
        """
        #!/usr/bin/env Rscript
        library(DESeq2)
        targets <- read.delim("$metadata",  header = TRUE)
        print(targets)
        directory = '/fsimb/groups-external/sfbcardosogr/ishita/maruthi_data/filtered_counts'
        group <- factor(targets\$group)
        donor <- factor(targets\$donor)
        ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = targets,
                                       directory = directory,
                                       design = ~group)
        dds <- DESeq(ddsHTSeq)
        res <- results(dds, alpha=0.01)
        write.table(res, "deg.tsv", sep="\t", row.names=FALSE, quote=FALSE)
        """
}

