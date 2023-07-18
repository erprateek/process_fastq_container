#!/bin/bash

set -euxo pipefail

TRIM_GALORE_OUTPUTS=~/trim_galore/
FASTQC_OUTPUTS_DIR=~/fastqc_outputs_post_trimming/

print_usage() {
  echo "Usage: $0 -1 <read1.fq.gz> -2 <read2.fq.gz>"
}



# Default values
read1=""
read2=""

# Parse arguments using getopt
options=$(getopt -o 1:2: -n "$0" -- "$@")
if [ $? -ne 0 ]; then
  print_usage
  exit 1
fi

eval set -- "$options"

while true; do
  case "$1" in
    -1)
      read1="$2"
      shift 2 ;;
    -2)
      read2="$2"
      shift 2 ;;
    --)
      shift
      break ;;
    *)
      print_usage
      exit 1 ;;
  esac
done

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
trim_galore --output $TRIM_GALORE_OUTPUTS --fastqc --fastqc_args "--outdir $FASTQC_OUTPUTS_DIR --extract" --paired $read1 $read2 
