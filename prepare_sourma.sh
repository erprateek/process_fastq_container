#!/bin/bash


mkdir -p /sourmash;
cd /sourmash;
curl -JLO https://osf.io/3a6gn/download;
curl -JLO https://osf.io/v3zmg/download;
sourmash sig summarize gtdb-rs207.genomic-reps.dna.k31.zip;
gunzip gtdb-rs207.taxonomy.csv.gz;
sourmash tax prepare -t gtdb-rs207.taxonomy.csv -o gtdb-rs207.taxonomy.sqldb -F sql
export SOURMASH_SQLDB=$(pwd)/gtdb-rs207.taxonomy.sqldb;
cd /;
