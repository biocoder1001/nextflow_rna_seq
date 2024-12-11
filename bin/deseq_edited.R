#!/usr/bin/env Rscript

library(org.Hs.eg.db)
library(optparse)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

alpha = 0.01
option_list <- list(
	make_option(c("-m", "--metadata"), type="character", default=NULL, metavar="path", help="metadata file."),
	make_option(c("-d", "--directory"), type="character", default=NULL, metavar="path", help="directory for the count files for the samples"),
	make_option(c("-r", "--ref"), type="character", default=NULL, metavar="path", help="reference group for the analysis"),
	make_option(c("-f", "--group_1"), type="character", default=NA, metavar="val", help="group_1 for the analysis"),
	make_option(c("-s", "--group_2"), type="character", default=NA, metavar="val", help="group_2 for the analysis"),
	make_option(c("-u", "--file_up_genes"), type="character", default="upregulated_genes.tsv", metavar="val", help="file name for upregulated genes for the analysis"),
	make_option(c("-l", "--file_down_genes"), type="character", default="downregulated_genes.tsv", metavar="val", help="file name for downregulated genes for the analysis")
)


opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$metadata)){
    print_help(opt_parser)
    stop("Please provide a metadata file.", call.=FALSE)
}


targets <- read.delim(opt$metadata,  header = TRUE)
directory = file.path(opt$directory, "filtered_counts")
if ("donor" %in% colnames(targets)) {
    targets$donor <- factor(targets$donor)
    group <- factor(targets$group)
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = targets,
                                       directory = directory,
                                       design = ~donor+group)
} else {
    targets$group <- factor(targets$group)
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = targets,
                                       directory = directory,
                                       design = ~group)
}


ddsHTSeq$group <- relevel(ddsHTSeq$group, ref=opt$ref)
dds <- DESeq(ddsHTSeq)
vsd <- vst(dds, blind=FALSE)
##quantification <- apply(fpm(dds), 2, function(x, y) 1e3 * x / y, gene.lengths[rownames(fpm(dds))])
##colnames(quantification) <- paste0(colnames(quantification),".robustFPKM")



##write.table(quantification, file=quantification.csv)
if (opt$group_2 == "NA") {
    res <- results(dds, contrast = list(opt$group_1))
} else {
    res <- results(dds, contrast = list(opt$group_1, opt$group_2))
}
print(colnames(res))

upregulated_genes <- subset(res, res$padj<alpha&res$log2FoldChange>0)
upregulated_genes["gene"] <- row.names(upregulated_genes)
downregulated_genes <- subset(res,  res$padj<alpha&res$log2FoldChange<0)
downregulated_genes$gene <- row.names(downregulated_genes)
upregulated_genes$gene_name <- mapIds(org.Hs.eg.db, keys=as.character(upregulated_genes$gene), keytype = 'ENSEMBL',column='SYMBOL', multiVals = "first")
downregulated_genes$gene_name <- mapIds(org.Hs.eg.db, keys=as.character(downregulated_genes$gene), keytype = 'ENSEMBL',column='SYMBOL', multiVals = "first")
write.table(downregulated_genes, file=opt$file_up_genes)
write.table(upregulated_genes, file=opt$file_down_genes)
pdf(file="DE_DESeq2.pdf")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$group
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
print(pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors))
print(plotPCA(vsd, intgroup = "group"))
dev.off()


