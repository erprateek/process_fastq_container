# process_fastq_container

This repository contains a solution for the Technical Coding Exercise, which involves processing paired-end sequencing data from an Illumina MiSeq library. The goal is to create a docker image (or other container type) that takes two input files in the fastq.gz format and generates two output files. The sequenced organism is a microbial isolate.

## Input Files
The input files provided for this exercise are two fastq.gz files representing paired-end sequencing data. These files contain the raw sequencing data from the Illumina MiSeq instrument.

## Output Files
The docker image is designed to generate the following output files:

1. **report.csv**: This is a comma-separated values (CSV) or text (TXT) report file that includes the following information:
   - Number of molecules sequenced
   - Sequencing configuration (read 1, index 1, index 2, and read 2 lengths)
   - Library type (Nextera or TruSeq)
   - Mean and standard deviation of library insert sizes after adapter removal
   - RefSeq or Genbank accession numbers of the sequenced organism

2. **plot.png**: This is a plot showing the instrument-reported Phred scores and alignment-inferred quality scores across all positions of read 1 and read 2. The x-axis represents positions across the read in units of base pairs, and the y-axis represents quality scores calculated as -10*log10(error rate). The plot can show Read 1 and Read 2 combined or as separate plots. It includes two lines representing the mean Phred qualities reported by the instrument and quality scores inferred from alignment error rates.

## Submission Instructions
Please follow the instructions below to build and run the Dockerfile, allocate sufficient memory, and estimate the execution time:

1. Clone this repository to your local machine.
2. Navigate to the repository's directory.
3. Build the Docker image using the provided Dockerfile and associated Python/Bash scripts. Make sure to create any necessary environment inside the Docker image.
4. Allocate enough memory to ensure smooth execution of the Docker image. The recommended memory allocation is [specify memory here].
5. Run the Docker image with the appropriate command, specifying the input files as arguments.
6. Wait for the execution to complete. The estimated execution time is [specify time here].
7. Retrieve the output files from the Docker image. The files you should find are **sample.bam** and either **report.txt** or **report.csv**.

Please make sure to zip the entire folder containing the Dockerfile, scripts, and output files. Submit the zipped folder by July 17, 2023.

## Contact
If you have any questions or need further clarification, please don't hesitate to reach out.
