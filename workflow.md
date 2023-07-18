# Order of steps to execute

## Installations



 ## Steps
 - Install miniconda
	 - Create environment tax_align: `conda create -n tax_align python=3.10`
 - install fastqc via `conda install -c conda-forge -c bioconda fastqc`
 - install cutadapt via `conda install -c conda-forge -c bioconda cutadapt`
 - Prepare sourmash:
	 - `conda install -c conda-forge -c bioconda sourmash parallel`
	 - `curl -JLO https://osf.io/3a6gn/download`
	 - `sourmash sig summarize gtdb-rs207.genomic-reps.dna.k31.zip`
	 - `curl -JLO https://osf.io/v3zmg/download`
	 - `gunzip gtdb-rs207.taxonomy.csv.gz`
	 - `sourmash tax prepare -t gtdb-rs207.taxonomy.csv     -o gtdb-rs207.taxonomy.sqldb -F sql`
 - Run sourmash on our reads
	 - `sourmash sketch dna -p k=31,abund /root/git/process_fastq_container/fastqs/read*.fastq -o sample.sig.gz --name unknown_sample`
	 - `sourmash gather sample.sig.gz gtdb-rs207.genomic-reps.dna.k31.zip --save-matches matches.zip`
 - For the organism identified in sourmash outputs, execute the cURL command to download reference:
	 - ` curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/<ORGANISM_ID>/download?include_annotation_type=GENOME_FASTA&filename=ORGANISM_ID.zip" -H "Accept: application/zip"`
	 - `unzip ORGANISM_ID.zip`
 - Install bwa-mem2:
	 - `git clone https://github.com/lh3/bwa.git`
	 - `cd bwa`
	 - `make`
 - Create reference index:
	 - `./bwa index ~/data/ncbi_dataset/data/GCF_000155815.1/GCF_000155815.1_ASM15581v1_genomic.fna`
 - Run Alignment:
	 - `./bwa mem ~/data/ncbi_dataset/data/GCF_000155815.1/GCF_000155815.1_ASM15581v1_genomic.fna ~/software/TrimGalore-0.6.10/read1_trimmed.fq ~/software/TrimGalore-0.6.10/read2_trimmed.fq > sample.sam`
 - Create alignment bam:
	 - `samtools view -bS sample.sam > sample.bam`

 
