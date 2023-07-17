FROM ubuntu:focal
LABEL maintainer="prateektandon@alumni.cmu.edu"
LABEL build_date="07/17/2023"

# install dependencies
RUN export DEBIAN_FRONTEND=noninteractive && \
apt-get update -y && \
apt-get install -y --no-install-recommends wget curl gunzip ca-certificates samtools libxml2 && \
apt-get clean && rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
bash -b -p ~/miniconda.sh

COPY environment.yml .
RUN conda env create -f environment.yml
SHELL ["conda", "run", "-n", "tax_align", "/bin/bash", "-c"]

# RUN conda install -c conda-forge -c bioconda fastqc cutadapt --yes
# RUN conda install -c conda-forge -c bioconda sourmash parallel && mkdir -p ~/sourmash && \

# Install bwa
# RUN git clone https://github.com/lh3/bwa.git ~/bwa && cd ~/bwa && make

# Prepare and install sourmash
cd ~/sourmash && curl -JLO https://osf.io/3a6gn/download && curl -JLO https://osf.io/v3zmg/download && \
sourmash sig summarize gtdb-rs207.genomic-reps.dna.k31.zip && \
gunzip gtdb-rs207.taxonomy.csv.gz && \
sourmash tax prepare -t gtdb-rs207.taxonomy.csv -o gtdb-rs207.taxonomy.sqldb -F sql

## Process data
RUN 