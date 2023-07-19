#!/bin/python

import os
import sys
import glob
import argparse
import gzip
import pandas as pd
import numpy as np

def get_organism_from_sourmash(sourmash_csv):
    df = pd.read_csv(sourmash_csv)
    name = df['name'].iloc[0]
    name_split = df['name'].split("=")[0]
    ncbi_id = name_split.split(" ")[0]
    organism_name = name_split.replace(ncbi_id, "").replace("strain","").lstrip().rstrip()
    return (ncbi_id, organism_name)

def get_isize_from_fq(fq):
    seq_index_base = 1
    seq_index_increment = 4
    isize_lens = []
    with gzip.open(fq, 'rb') as fh:
        fqlines = fh.readlines()
        for i in range(0, len(fqlines), seq_index_increment):
            line_to_process = fqlines[i].strip()
            isize_lens.append(len(line_to_process))
    return isize_lens
        
def get_insert_size_mean_and_std(fq1, fq2):
    isize_mean = 0
    isize_std = 0
    isize_dist_fq1 = get_isize_from_fq(fq1)
    isize_dist_fq2 = get_isize_from_fq(fq2)
    mean_fq1 = np.mean(isize_dist_fq1)
    mean_fq2 = np.mean(isize_dist_fq2)
    std_fq1 = np.std(isize_dist_fq1)
    std_fq2 = np.std(isize_dist_fq2)
    mean_isize = np.mean(isize_dist_fq1+isize_dist_fq2)
    std_isize = np.std(isize_dist_fq1+isize_dist_fq2)
    return (isize_mean, isize_std)

def get_read_len_and_ilen_for_fq(fq):
    seq_len = 0
    ilen = 0
    with gzip.open(fq, 'rb') as fh:
        fqlines = fh.readlines()
        header = fqlines[0].strip()
        sequence = fqlines[1].strip()
        index_seq = header.split(":")[-1]
        seq_len = len(sequence)
        ilen = len(index_seq)
    return (seq_len, ilen)

def parse_trimgalore_outputs(trimgalore_output_dir):
    reports = sorted(glob.glob(os.path.join(trimgalore_output_dir, "report.txt")))
    content = []
    with open(reports[0]) as tfile:
        content = tfile.readlines()
    adapter_sequence = ""
    library_type = "Unknown"
    for line in content:
        if line.startswith("Adapter sequence"):
            adapter_sequence = line.split(" ")[2].replace("'","")
            if "nextera" in line.lower():
                library_type = "Nextera"
            elif "truseq" in line.lower():
                library_type = "Truseq"
    return (adapter_sequence, library_type)

def get_sequencing_configuration(fq1, fq2):
    r1_len, i1_len = get_read_len_and_ilen_for_fq(fq1)
    r2_len, i2_len = get_read_len_and_ilen_for_fq(fq2)
    return (r1_len, r2_len, i1_len, i2_len)

def count_reads_in_fastq(fastq_file):
    count = 0
    with gzip.open(fastq_file, 'r') as file:
        for line in file:
            if line.startswith('@'):
                count += 1
    return count

def get_num_molecules(read1_fastq, read2_fastq):
    count_read1 = count_reads_in_fastq(read1_fastq)
    count_read2 = count_reads_in_fastq(read2_fastq)
    return count_read1+count_read2
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",'--fastqc_outputs_dir')
    parser.add_argument("-1","--read1")
    parser.add_argument("-2","--read2")
    parser.add_argument('-3','--read1_trimmed')
    parser.add_argument('-4','--read2_trimmed')
    parser.add_argument("-t"."--trimgalore_output_dir")
    parser.add_argument('-s','--sourmash_csv')
    #parser.add_argument("-b",'--sample_bam')
    parser.add_argument("-o","--output_file")  
    args = parser.parse_args()
    
    isize_mean, isize_std = get_insert_size_mean_and_std(args.read1_trimmed, args.read2_trimmed)
    ncbi_id, name = get_organism_from_sourmash(args.sourmash_csv)
    r1_len, r2_len, i1_len, i2_len = get_sequencing_configuration(args.read1, args.read2)
    adapter, library = parse_trimgalore_outputs(args.trimgalore_output_dir)
    nmolecules = get_num_molecules(args.read1, args.read2)
    row = [nmolecules, r1_len, r2_len, i1_len, i2_len, library, isize_mean, isize_std, ncbi_id]
    cols = [
        'Num_molecules',
        'R1_length',
        'R2_length',
        'Index1_length',
        'Index2_length',
        'Library',
        'Insert_size_mean',
        'Insert_size_standard_deviation',
        'NCBI_ID'
    ]
    df = pd.DataFrame(row, columns=cols)
    df.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()