#!/bin/bash

set -euxo pipefail

TRIM_GALORE_OUTPUTS=/trim_galore/
FASTQC_OUTPUTS_DIR=/fastqc_outputs_post_trimming/
SOURMASH_OUTPUTS=/sourmash_outputs/
REFERENCES_DIR=/references/
PIPELINE_OUTPUTS=/outputs
print_usage() {
  echo "Usage: $0 <read1.fq.gz> <read2.fq.gz>"
}

# Check if both read1 and read2 were provided
if [ -z "$read1" ] || [ -z "$read2" ]; then
  echo "Error: Both read1 and read2 files are required."
  print_usage
  exit 1
fi

# Check if the provided files exist
if [ ! -f "$read1" ] || [ ! -f "$read2" ]; then
  echo "Error: One or both of the input files do not exist."
  exit 1
fi

# Add your desired actions here, e.g., perform processing on read1 and read2

# Example: Print the file names
echo "Obtained read1 file: $read1"
echo "Obtained read2 file: $read2"
mkdir -p $FASTQC_OUTPUTS_DIR
mkdir -p $TRIM_GALORE_OUTPUTS
mkdir -p $SOURMASH_OUTPUTS
mkdir -p $PIPELINE_OUTPUTS

trim_galore --output $TRIM_GALORE_OUTPUTS --fastqc --fastqc_args "--outdir $FASTQC_OUTPUTS_DIR --extract" --paired $read1 $read2
R1_TRIMMED=$TRIM_GALORE_OUTPUTS/

SOURMASH_SIGNATURE_FILE=sample.sig.gz
SOURMASH_K=31
SOURMASH_OUTPUT_CSV=sourmash_output.csv
sourmash sketch dna -p k=$SOURMASH_K,abund $TRIM_GALORE_OUTPUTS/read*.fastq.gz -o $SOURMASH_OUTPUTS/$SOURMASH_SIGNATURE_FILE --name new_sample
sourmash gather $SOURMASH_SIGNATURE_FILE gtdb-rs207.genomic-reps.dna.k$SOURMASH_K.zip --output $SOURMASH_OUTPUTS/$SOURMASH_OUTPUT_CSV


## Download reference
mkdir -p $REFERENCES_DIR
NCBI_ID=$(cut -d',' -f 10 $SOURMASH_OUTPUT_CSV | grep -v name | head -1 | sed 's/"//g' | cut -d " " -f 1);
curl -o $REFERENCES_DIR -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$NCBI_ID/download?include_annotation_type=GENOME_FASTA&filename=${NCBI_ID}.zip" -H "Accept: application/zip" 
unzip $REFERENCES_DIR/$NCBI_ID.zip
REFERENCE=$REFERENCES_DIR/ncbi_dataset/data/$NCBI_ID/$NCBI_ID_genomic.fna

## BWA Index
bwa index $REFERENCE

## Run Alignment
bwa mem $REFERENCE $TRIM_GALORE_OUTPUTS/*1.fq.gz $TRIM_GALORE_OUTPUTS/*2.fq.gz | samtools view -bS $PIPELINE_OUTPUTS/sample.bam

## Create plot
python generate_report.py -f $FASTQC_OUTPUTS_DIR -t $TRIM_GALORE_OUTPUTS -s $SOURMASH_OUTPUTS -b $PIPELINE_OUTPUTS/sample.bam -o $PIPELINE_OUTPUTS/report.csv
#python create_plot.py -b $PIPELINE_OUTPUTS/sample.bam -1 $read1 -2 $read2 -o $PIPELINE_OUTPUTS/plot.png
