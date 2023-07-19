# Use the miniconda3 base image as the starting point.
# This base image includes the Miniconda package manager, which is useful for managing Python packages and environments.
FROM continuumio/miniconda3

# Set the maintainer's email and build date as labels for the image.
LABEL maintainer="prateektandon@alumni.cmu.edu"
LABEL build_date="07/19/2023"

# Install system dependencies using apt-get.
RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update -y && \
    apt-get install -y --no-install-recommends wget curl ca-certificates samtools libxml2 unzip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy the 'environment.yml' file into the image.
# This file contains the Conda environment specification.
COPY environment.yml .

# Create a Conda environment named 'tax_align' and install the required packages.
RUN conda env create -f environment.yml

# Configure the 'bashrc' to automatically activate the 'tax_align' Conda environment when starting a new shell.
RUN echo "conda activate tax_align" >> ~/.bashrc

# Set the shell to use 'bash' with login options and running commands as a login shell.
SHELL ["/bin/bash", "--login", "-c"]

# Verify that the required Python package 'numpy' is installed in the 'tax_align' Conda environment.
RUN echo "Make sure modules are is installed:"
RUN python -c "import numpy"

# Copy the 'prepare_sourma.sh' script into the image.
COPY prepare_sourma.sh .

# Run the 'prepare_sourma.sh' script to prepare and install sourmash.
RUN bash prepare_sourma.sh

# Copy the 'execute_workflow.sh', 'generate_report.py', and 'create_plot.py' files into the image.
COPY execute_workflow.sh .
COPY generate_report.py .
COPY create_plot.py .

# Move the script files to '/usr/local/bin' and set them as executable.
RUN cp execute_workflow.sh /usr/local/bin/execute_workflow.sh && chmod a+x /usr/local/bin/execute_workflow.sh;
RUN cp generate_report.py /usr/local/bin/generate_report.py && chmod a+x /usr/local/bin/generate_report.py;
RUN cp create_plot.py /usr/local/bin/create_plot.py && chmod a+x /usr/local/bin/create_plot.py;
