#FROM ubuntu:focal
FROM continuumio/miniconda3
LABEL maintainer="prateektandon@alumni.cmu.edu"
LABEL build_date="07/17/2023"

# install dependencies
RUN export DEBIAN_FRONTEND=noninteractive && \
apt-get update -y && \
apt-get install -y --no-install-recommends wget curl ca-certificates samtools libxml2 && \
apt-get clean && rm -rf /var/lib/apt/lists/*

#RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#bash ~/miniconda.sh -b -p $HOME/miniconda

COPY environment.yml .
RUN conda env create -f environment.yml
RUN echo "conda activate tax_align" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]
#SHELL ["conda", "run", "-n", "tax_align", "/bin/bash", "-c"]

RUN echo "Make sure modules are is installed:"
RUN python -c "import numpy"

# Prepare and install sourmash
COPY prepare_sourma.sh .
RUN bash prepare_sourma.sh
