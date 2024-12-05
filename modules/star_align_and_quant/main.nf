process alignment_quant {
        publishDir  "${params.basedir}/aligned_bam_files", mode: 'copy'
        input:
        path reads
        path index
        path gtf

        output:
        path "bam_files/*.bam", emit: bam_files
        path "count_files/*.ReadsPerGene.out.tab", emit: count_matrix

        script:
        """
        mkdir -p bam_files
        mkdir -p count_files
        STAR --runThreadN 4 \
         --quantMode GeneCounts\
         --runMode alignReads\
         --outFilterMismatchNmax 999\
         --sjdbOverhang 49\
         --sjdbGTFfile $gtf\
         --genomeDir $index \
         --readFilesIn $reads \
         --readFilesCommand zcat \
         --outFileNamePrefix bam_files/${reads[0].baseName}. \
         --outSAMtype BAM SortedByCoordinate
         mv bam_files/*.ReadsPerGene.out.tab count_files/
        """
}

