#!/bin/bash
source ~/.bashrc
set -euxo pipefail

TRIM_GALORE_OUTPUTS=/trim_galore/
FASTQC_OUTPUTS_DIR=/fastqc_outputs_post_trimming/
SOURMASH_OUTPUTS=/sourmash_outputs/
REFERENCES_DIR=/references/
PIPELINE_OUTPUTS=/outputs

# Check if the correct number of parameters is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <read1.fastq> <read2.fastq>"
    exit 1
fi

# Check if the provided files exist
if [ ! -f "$1" ]; then
    echo "Error: File '$1' not found."
    exit 1
fi

if [ ! -f "$2" ]; then
    echo "Error: File '$2' not found."
    exit 1
fi

read1=$1
read2=$2
# Add your desired actions here, e.g., perform processing on read1 and read2

# Example: Print the file names
echo "Obtained read1 file: $read1"
echo "Obtained read2 file: $read2"
mkdir -p $FASTQC_OUTPUTS_DIR
mkdir -p $TRIM_GALORE_OUTPUTS
mkdir -p $SOURMASH_OUTPUTS
mkdir -p $PIPELINE_OUTPUTS

trim_galore --output $TRIM_GALORE_OUTPUTS --fastqc --fastqc_args "--outdir $FASTQC_OUTPUTS_DIR --extract" --paired $read1 $read2
R1_TRIMMED=$TRIM_GALORE_OUTPUTS/*1.*.gz
R2_TRIMMED=$TRIM_GALORE_OUTPUTS/*2.*.gz

SOURMASH_SIGNATURE_FILE=sample.sig.gz
SOURMASH_K=31
SOURMASH_OUTPUT_CSV=sourmash_output.csv
sourmash sketch dna -p k=$SOURMASH_K,abund $R1_TRIMMED $R2_TRIMMED -o $SOURMASH_OUTPUTS/$SOURMASH_SIGNATURE_FILE --name new_sample
sourmash gather $SOURMASH_OUTPUTS/$SOURMASH_SIGNATURE_FILE /sourmash/gtdb-rs207.genomic-reps.dna.k$SOURMASH_K.zip --output $SOURMASH_OUTPUTS/$SOURMASH_OUTPUT_CSV

## Download reference
mkdir -p $REFERENCES_DIR
NCBI_ID=$(cut -d',' -f 10 $SOURMASH_OUTPUTS/$SOURMASH_OUTPUT_CSV | grep -v name | head -1 | sed 's/"//g' | cut -d " " -f 1);
curl -o $REFERENCES_DIR/$NCBI_ID.zip -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$NCBI_ID/download?include_annotation_type=GENOME_FASTA&filename=${NCBI_ID}.zip" -H "Accept: application/zip"
cd $REFERENCES_DIR
unzip $NCBI_ID.zip
REFERENCE=$REFERENCES_DIR/ncbi_dataset/data/${NCBI_ID}/${NCBI_ID}*_genomic.fna
cd -
## BWA Index
bwa index $REFERENCE

## Run Alignment
bwa mem $REFERENCE $TRIM_GALORE_OUTPUTS/*1.fq.gz $TRIM_GALORE_OUTPUTS/*2.fq.gz | samtools view -bS - > $PIPELINE_OUTPUTS/sample.bam

## Create plot
#python generate_report.py -f $FASTQC_OUTPUTS_DIR -t $TRIM_GALORE_OUTPUTS -s $SOURMASH_OUTPUTS -b $PIPELINE_OUTPUTS/sample.bam -o $PIPELINE_OUTPUTS/report.csv
python generate_report.py -1 $read1 -2 $read2 -3 $TRIM_GALORE_OUTPUTS/*1.fq.gz -4 $TRIM_GALORE_OUTPUTS/*2.fq.gz -t $TRIM_GALORE_OUTPUTS -s $SOURMASH_OUTPUTS/$SOURMASH_OUTPUT_CSV -o $PIPELINE_OUTPUTS/report.csv
#python generate_report.py -1 fastqs/read1.fastq.gz -2 fastqs/read2.fastq.gz -3 /root/trimgalore_outputs/read1_val_1.fq.gz -4 /root/trimgalore_outputs/read2_val_2.fq.gz -t /root/trimgalore_outputs/ -s /root/software/sourmash/sourmash_gather_output.csv -o report.csv
samtools view $PIPELINE_OUTPUTS/sample.bam | cut -f 11 > /bam_qc_output.txt
python create_plot.py -b /bam_qc_output.txt -1 $read1 -2 $read2 -o $PIPELINE_OUTPUTS/plot.png
