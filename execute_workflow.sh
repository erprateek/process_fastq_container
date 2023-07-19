#!/bin/bash

# Load user-defined settings from the bashrc file.
source ~/.bashrc

# Enable strict error handling for the script.
set -euxo pipefail

# Set up directories for various outputs and references.
TRIM_GALORE_OUTPUTS=/trim_galore/
FASTQC_OUTPUTS_DIR=/fastqc_outputs_post_trimming/
SOURMASH_OUTPUTS=/sourmash_outputs/
REFERENCES_DIR=/references/
PIPELINE_OUTPUTS=/outputs

# Check if the correct number of parameters is provided.
if [ $# -ne 2 ]; then
    echo "Usage: $0 <read1.fastq> <read2.fastq>"
    exit 1
fi

# Check if the provided files exist.
if [ ! -f "$1" ]; then
    echo "Error: File '$1' not found."
    exit 1
fi

if [ ! -f "$2" ]; then
    echo "Error: File '$2' not found."
    exit 1
fi

# Assign read1 and read2 file paths based on the provided parameters.
read1=$1
read2=$2

# Create necessary directories for the pipeline outputs.
mkdir -p $FASTQC_OUTPUTS_DIR
mkdir -p $TRIM_GALORE_OUTPUTS
mkdir -p $SOURMASH_OUTPUTS
mkdir -p $PIPELINE_OUTPUTS

# Run 'trim_galore' to perform trimming on read1 and read2, and generate FastQC reports.
trim_galore --output $TRIM_GALORE_OUTPUTS --fastqc --fastqc_args "--outdir $FASTQC_OUTPUTS_DIR --extract" --paired $read1 $read2
R1_TRIMMED=$TRIM_GALORE_OUTPUTS/*1.*.gz
R2_TRIMMED=$TRIM_GALORE_OUTPUTS/*2.*.gz

# Prepare and run 'sourmash' to create a signature file and gather taxonomic information.
SOURMASH_SIGNATURE_FILE=sample.sig.gz
SOURMASH_K=31
SOURMASH_OUTPUT_CSV=sourmash_output.csv
sourmash sketch dna -p k=$SOURMASH_K,abund $R1_TRIMMED $R2_TRIMMED -o $SOURMASH_OUTPUTS/$SOURMASH_SIGNATURE_FILE --name new_sample
sourmash gather $SOURMASH_OUTPUTS/$SOURMASH_SIGNATURE_FILE /sourmash/gtdb-rs207.genomic-reps.dna.k$SOURMASH_K.zip --output $SOURMASH_OUTPUTS/$SOURMASH_OUTPUT_CSV

# Download the reference genome data for taxonomic alignment.
mkdir -p $REFERENCES_DIR
NCBI_ID=$(cut -d',' -f 10 $SOURMASH_OUTPUTS/$SOURMASH_OUTPUT_CSV | grep -v name | head -1 | sed 's/"//g' | cut -d " " -f 1);
curl -o $REFERENCES_DIR/$NCBI_ID.zip -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$NCBI_ID/download?include_annotation_type=GENOME_FASTA&filename=${NCBI_ID}.zip" -H "Accept: application/zip"
cd $REFERENCES_DIR
unzip $NCBI_ID.zip
REFERENCE=$REFERENCES_DIR/ncbi_dataset/data/${NCBI_ID}/${NCBI_ID}*_genomic.fna
cd -

# Create BWA index for the downloaded reference genome.
bwa index $REFERENCE

# Perform read alignment using BWA-MEM and generate a BAM file.
bwa mem $REFERENCE $TRIM_GALORE_OUTPUTS/*1.fq.gz $TRIM_GALORE_OUTPUTS/*2.fq.gz | samtools view -bS - > $PIPELINE_OUTPUTS/sample.bam

# Generate a report using Python script 'generate_report.py' based on the processed data.
python generate_report.py -1 $read1 -2 $read2 -3 $TRIM_GALORE_OUTPUTS/*1.fq.gz -4 $TRIM_GALORE_OUTPUTS/*2.fq.gz -t $TRIM_GALORE_OUTPUTS -s $SOURMASH_OUTPUTS/$SOURMASH_OUTPUT_CSV -o $PIPELINE_OUTPUTS/report.csv

# Extract the mapping quality (MAPQ) from the BAM file and generate a plot using 'create_plot.py'.
samtools view $PIPELINE_OUTPUTS/sample.bam | cut -f 11 > /bam_qc_output.txt
python create_plot.py -b /bam_qc_output.txt -1 $TRIM_GALORE_OUTPUTS/*1.fq.gz -2 $TRIM_GALORE_OUTPUTS/*2.fq.gz -o $PIPELINE_OUTPUTS/plot.png
