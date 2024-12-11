#!/usr/bin/env nextflow

/*
* This is a rnaseq pipeline
*/

nextflow.enable.dsl=2

log.info """\
        RNASEQ-NF-PIPELINE
        ==================
        ISHITA | CARDOSO LAB | TECHNICAL UNIVERSITY OF DARMSTADT
        """


def helpMessage(){
        log.info """
                Usage:
                nextflow run main-nf --index --reads --gtf_file --samples --basedir
                Mandatory arguments:
                --index:path to index file for star
                --reads:path to reads (pipeline is for single end
                --gtf_file:path to gtf file for the reference genome
                --samples:sample file with sample names that are to be used for differential expression analysis
                --basedir:path to the folder where the script is located
                """

}

if (params.help) {
    helpMessage()
    exit 0
}

// Define channels
pairedEndRegex = params.inputDir + "/*_{1,2}.fastq"
SERegex = params.inputDir + "/*[!12].fastq"
pairFiles = Channel.fromFilePairs(pairedEndRegex)
singleFiles = Channel.fromFilePairs(SERegex, size: 1){ file -> file.baseName.replaceAll(/.fastq/,"")}
singleFiles.mix(pairFiles)
.set { reads_ch }
//reads_ch = Channel.fromFilePairs(params.reads,checkIfExists: true, size:params.singleEnd ? 1:2)
index_ch = Channel.value(params.index)
gtf_ch = Channel.value(params.gtf_file)
sample_selection = Channel.value(params.samples)
ref = Channel.value(params.ref)



// import modules
include { fastqc } from './modules/fastqc'
include { cutadapt } from './modules/cutadapt'
include { alignment_quant } from './modules/star_align_and_quant'
include { counts_mat_processing } from './modules/count_processing'
include { sample_file_extract } from './modules/sample_file_extract_for_deseq'
include { metadata_prep} from './modules/metadata_gen'
include { deseq_analysis} from './modules/deseq_analysis'


workflow {

  fastqc(reads_ch)
  cutadapt(reads_ch)
  alignment_quant(cutadapt.out,index_ch,gtf_ch)
  counts_mat_processing(alignment_quant.out.count_matrix)
  counts=sample_file_extract(counts_mat_processing.out,sample_selection).filter{file ->file.name != 'empty.tsv'}
  counts.view()
  count_files_list = counts.collect { file -> file.name }
  metadata_prep(count_files_list,sample_selection)
  deseq_analysis(metadata_prep.out, ref, params.group_1, params.group_2)
}
