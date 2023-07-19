#FROM ubuntu:focal
FROM continuumio/miniconda3
LABEL maintainer="prateektandon@alumni.cmu.edu"
LABEL build_date="07/19/2023"

# install dependencies
RUN export DEBIAN_FRONTEND=noninteractive && \
apt-get update -y && \
apt-get install -y --no-install-recommends wget curl ca-certificates samtools libxml2 unzip && \
apt-get clean && rm -rf /var/lib/apt/lists/*

COPY environment.yml .
RUN conda env create -f environment.yml
RUN echo "conda activate tax_align" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

RUN echo "Make sure modules are is installed:"
RUN python -c "import numpy"

# Prepare and install sourmash
COPY prepare_sourma.sh .
RUN bash prepare_sourma.sh

COPY execute_workflow.sh .
COPY generate_report.py .
COPY create_plot.py .

RUN cp execute_workflow.sh /usr/local/bin/execute_workflow.sh && chmod a+x /usr/local/bin/execute_workflow.sh;
RUN cp generate_report.py /usr/local/bin/generate_report.py && chmod a+x /usr/local/bin/generate_report.py;
RUN cp create_plot.py /usr/local/bin/create_plot.py && chmod a+x /usr/local/bin/create_plot.py;
