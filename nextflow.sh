#!/bin/bash 

##############################
#       Job blueprint        #
##############################

# Give your job a name, so you can recognize it in the queue overview
#SBATCH --job-name=rnaseq

# Define, how many nodes you need. Here, we ask for 1 node.
# Each node has 16 or 20 CPU cores.
#SBATCH --nodes=2
#--ntasks=4
# You can further define the number of tasks with --ntasks-per-*
# See "man sbatch" for details. e.g. --ntasks=4 will ask for 4 cpus.

# Define, how long the job will run in real time. This is a hard cap meaning
# that if the job runs longer than what is written here, it will be
# force-stopped by the server. If you make the expected time too long, it will
# take longer for the job to start. Here, we say the job will take 5 minutes.
#              d-hh:mm:ss
#SBATCH --time=1-12:30:00

# Define the partition on which the job shall run. May be omitted.
#SBATCH --partition long

# How much memory you need.
# --mem will define memory per node and
# --mem-per-cpu will define memory per CPU/core. Choose one of those.
##SBATCH --mem-per-cpu=1500MB
#SBATCH --mem=900GB    # this one is not in effect, due to the double hash

# Turn on mail notification. There are many possible self-explaining values:
# NONE, BEGIN, END, FAIL, ALL (including all aforementioned)
# For more values, check "man sbatch"
#SBATCH --mail-type=END,FAIL
##module load nextflow/23.10.0
module load cutadapt/4.0
module load R/Bioconductor_3.18_singularity 
module load fastqc/0.11.9 
module load star/2.7.3a 
module load nextflow/23.10.0
module load jdk/18.0.2.1 
nextflow run nextflow_workflow.nf -c nextflow.config --reads "/fsimb/groups-external/sfbcardosogr/ishita/nextflow_pipelines/fastq_files/*.fastq.gz"  --index /fsimb/groups-external/sfbcardosogr/ishita/ensembl/grch38/canonical/index/star/2.7.3a  --gtf_file /fsimb/groups-external/sfbcardosogr/ishita/ensembl/grch38/canonical/annotation/Homo_sapiens.GRCh38.98.gtf  --basedir /fsimb/groups-external/sfbcardosogr/ishita/nextflow-pipelines  -with-trace -with-report -with-timeline
