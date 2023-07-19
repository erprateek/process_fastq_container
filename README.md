# Taxonomic Alignment Docker Image

This Docker image contains a pre-configured environment for running a taxonomic alignment workflow. The image is built on top of the Miniconda3 base image and includes all the necessary dependencies to execute the workflow.

## Prerequisites

Before you proceed, make sure you have the following software installed on your system:

1. Docker: Install Docker from the official website based on your operating system. Refer to [Docker Installation Guide](https://docs.docker.com/get-docker/) for detailed instructions.

## Building the Docker Image

To build the Docker image from the provided Dockerfile, follow these steps:

1. Clone the repository containing the Dockerfile and associated scripts.

2. Open a terminal or command prompt and navigate to the root directory of the cloned repository.

3. Run the following Docker command to build the image:

```bash
docker build -t tax_align_image .
```
## Running the docker container with example fastqs:
```bash
docker run --rm -v <PATH_TO_FASTQ>/read1.fastq.gz:/read1.fastq.gz \
-v <PATH_TO_FASTQ>/read2.fastq.gz:/read2.fastq.gz \
-v <PATH_TO_OUTPUTS_DIR_ON_HOST>:/outputs tax_align_image bash execute_workflow.sh read1.fastq.gz read2.fastq.gz
```