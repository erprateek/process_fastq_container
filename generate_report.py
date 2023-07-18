#!/bin/python

import os
import sys
import glob
import argparse
import pandas as pd

def get_organism_from_sourmash(sourmash_csv):
    df = pd.read_csv(sourmash_csv)
    name = df['name'].iloc[0]
    name_split = df['name'].split("=")[0]
    ncbi_id = name_split.split(" ")[0]
    organism_name = name_split.replace(ncbi_id, "").replace("strain","").lstrip().rstrip()
    return (ncbi_id, organism_name)

def get_insert_size_mean_and_std():
    isize_mean = 0
    isize_std = 0
    return (isize_mean, isize_std)

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

def get_sequencing_configuration():
    r1_len = 0
    r2_len = 0
    i1_len = 0
    i2_len = 0
    return (r1_len, r2_len, i1_len, i2_len)

def get_num_molecules():
    num_molecules = 0
    return num_molecules    
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",'--fastqc_outputs_dir')
    parser.add_argument("-t"."--trimgalore_output_dir")
    parser.add_argument('-s','--sourmash_csv')
    parser.add_argument("-b",'--sample_bam')
    parser.add_argument("-o","--output_file")
    
    args = parser.parse_args()
    
    isize_mean, isize_std = get_insert_size_mean_and_std()
    ncbi_id, name = get_organism_from_sourmash(args.sourmash_csv)
    r1_len, r2_len, i1_len, i2_len = get_sequencing_configuration()
    adapter, library = parse_trimgalore_outputs(args.trimgalore_output_dir)
    nmolecules = get_num_molecules()
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