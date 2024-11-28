
#!/usr/bin/env Rscript
args = commandArgs(T)
alpha = 0.01
library(DESeq2)
library(org.Hs.eg.db)

targets <- read.delim(args[1],  header = TRUE)
directory = '/fsimb/groups-external/sfbcardosogr/ishita/nextflow-pipelines/filtered_counts'
group <- factor(targets$group)
donor <- factor(targets$donor)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = targets,
                                       directory = directory,
                                       design = ~group)
dds <- DESeq(ddsHTSeq)
res <- results(dds)
res1 <- results(dds)
upregulated_genes <- subset(res1, res1$padj<alpha&res1$log2FoldChange>0)
upregulated_genes["gene"] <- row.names(upregulated_genes)
downregulated_genes <- subset(res1,  res1$padj<alpha&res1$log2FoldChange<0)
downregulated_genes$gene <- row.names(downregulated_genes)
##upregulated_genes$gene_name <- mapIds(org.Hs.eg.db, keys=as.character(upregulated_genes$gene), keytype = 'ENSEMBL',column='SYMBOL', multiVals = "first")
##downregulated_genes$gene_name <- mapIds(org.Hs.eg.db, keys=as.character(downregulated_genes$gene), keytype = 'ENSEMBL',column='SYMBOL', multiVals = "first")
##df_list <- list(upregulated_genes, downregulated_genes)
write.table(downregualated_genes, file="up_and_down_regulated_genes_single_ended.tsv")
