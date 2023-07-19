# Import required Snakemake modules
from snakemake.io import expand

# Define the input and output files/directories
TRIM_GALORE_OUTPUTS = "/trim_galore/"
FASTQC_OUTPUTS_DIR = "/fastqc_outputs_post_trimming/"
SOURMASH_OUTPUTS = "/sourmash_outputs/"
REFERENCES_DIR = "/references/"
PIPELINE_OUTPUTS = "/outputs"

# Define the rule to check the input files
rule check_input_files:
    input:
        read1="{read1}",
        read2="{read2}"
    shell:
        """
        if [ ! -f "{input.read1}" ]; then
            echo "Error: File '{input.read1}' not found."
            exit 1
        fi

        if [ ! -f "{input.read2}" ]; then
            echo "Error: File '{input.read2}' not found."
            exit 1
        fi
        """

# Rule to perform read trimming using Trim Galore
rule trim_galore:
    input:
        read1="{read1}",
        read2="{read2}"
    output:
        trimmed_R1="{TRIM_GALORE_OUTPUTS}/{sample}_1.fq.gz",
        trimmed_R2="{TRIM_GALORE_OUTPUTS}/{sample}_2.fq.gz",
        fastqc_report="{FASTQC_OUTPUTS_DIR}/{sample}_fastqc.html"
    shell:
        """
        mkdir -p {FASTQC_OUTPUTS_DIR}
        mkdir -p {TRIM_GALORE_OUTPUTS}
        trim_galore --output {TRIM_GALORE_OUTPUTS} --fastqc --fastqc_args "--outdir {FASTQC_OUTPUTS_DIR} --extract" --paired {input.read1} {input.read2}
        """

# Rule to create Sourmash signature files
rule sourmash:
    input:
        R1="{TRIM_GALORE_OUTPUTS}/{sample}_1.fq.gz",
        R2="{TRIM_GALORE_OUTPUTS}/{sample}_2.fq.gz"
    output:
        signature="{SOURMASH_OUTPUTS}/{sample}.sig.gz"
    params:
        k=31
    shell:
        """
        sourmash sketch dna -p k={params.k},abund {input.R1} {input.R2} -o {output.signature} --name {wildcards.sample}
        """

# Rule to perform Sourmash gather
rule sourmash_gather:
    input:
        signature="{SOURMASH_OUTPUTS}/{sample}.sig.gz"
    output:
        csv="{SOURMASH_OUTPUTS}/{sample}_gather_output.csv"
    params:
        k=31
    shell:
        """
        sourmash gather {input.signature} /sourmash/gtdb-rs207.genomic-reps.dna.k{params.k}.zip --output {output.csv}
        """

# Rule to download reference data
rule download_reference:
    output:
        zip="{REFERENCES_DIR}/{ncbi_id}.zip"
    shell:
        """
        NCBI_ID=$(cut -d',' -f 10 {SOURMASH_OUTPUTS}/{wildcards.sample}_gather_output.csv | grep -v name | head -1 | sed 's/"//g' | cut -d " " -f 1)
        curl -o {output.zip} -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{wildcards.sample}/download?include_annotation_type=GENOME_FASTA&filename=${NCBI_ID}.zip" -H "Accept: application/zip"
        """

# Rule to unzip reference data
rule unzip_reference:
    input:
        zip="{REFERENCES_DIR}/{ncbi_id}.zip"
    output:
        fasta="{REFERENCES_DIR}/{ncbi_id}_genomic.fna"
    shell:
        """
        unzip -o {input.zip} -d {REFERENCES_DIR}
        cd {REFERENCES_DIR}/ncbi_dataset/data/{ncbi_id}
        mv *genomic.fna ../../{ncbi_id}_genomic.fna
        cd -
        """

# Rule to create BWA index
rule bwa_index:
    input:
        fasta="{REFERENCES_DIR}/{ncbi_id}_genomic.fna"
    output:
        bwa_index="{REFERENCES_DIR}/{ncbi_id}_genomic.fna.bwt"
    shell:
        """
        bwa index {input.fasta}
        """

# Rule to perform BWA alignment
rule bwa_mem:
    input:
        fasta="{REFERENCES_DIR}/{ncbi_id}_genomic.fna",
        R1="{TRIM_GALORE_OUTPUTS}/{sample}_1.fq.gz",
        R2="{TRIM_GALORE_OUTPUTS}/{sample}_2.fq.gz"
    output:
        bam="{PIPELINE_OUTPUTS}/{sample}.bam"
    shell:
        """
        bwa mem {input.fasta} {input.R1} {input.R2} | samtools view -bS - > {output.bam}
        """

# Rule to generate report
rule generate_report:
    input:
        read1="{read1}",
        read2="{read2}",
        R1="{TRIM_GALORE_OUTPUTS}/{sample}_1.fq.gz",
        R2="{TRIM_GALORE_OUTPUTS}/{sample}_2.fq.gz",
        csv="{SOURMASH_OUTPUTS}/{sample}_gather_output.csv"
    output:
        report="{PIPELINE_OUTPUTS}/{sample}_report.csv"
    shell:
        """
        python generate_report.py -1 {input.read1} -2 {input.read2} -3 {input.R1} -4 {input.R2} -t {TRIM_GALORE_OUTPUTS} -s {input.csv} -o {output.report}
        """

# Rule to create plot
rule create_plot:
    input:
        bam="{PIPELINE_OUTPUTS}/{sample}.bam"
    output:
        plot="{PIPELINE_OUTPUTS}/{sample}_plot.png"
    shell:
        """
        samtools view {input.bam} | cut -f 11 > {PIPELINE_OUTPUTS}/{sample}_bam_qc_output.txt
        python create_plot.py -b {PIPELINE_OUTPUTS}/{sample}_bam_qc_output.txt -1 {read1} -2 {read2} -o {output.plot}
        """
