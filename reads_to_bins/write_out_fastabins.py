#!/bin/bash 
import argparse
import sys
import os 
import re
from Bio import SeqIO
from collections import defaultdict

### Usage 

parser = argparse.ArgumentParser(description='''
   Takes a .clusters.tsv file from VAMB and subsets the fasta file with contigs used in VAMB
   ''')
parser.add_argument('--c', help='Cluster file')
parser.add_argument('--f', help='Fasta file of Encoded contigs')
parser.add_argument('--d', help='Out Directory of fasta files', type=str, default=None)


# parse the command line
args = parser.parse_args()

### Get Cluster - Contig associations 
CLUSTER_FILE=args.c
FASTA_FILE=args.f
OUTDIR=args.d

os.system('mkdir -p {}'.format(OUTDIR))

### Cluster file needs to contain
## cluster_1 c1
## cluster_1 c2
## etc.
print('Establishing Cluster and Contig relations')
contig_clusters = defaultdict(dict)
with open(CLUSTER_FILE,'r') as infile:
    for line in infile:
        line = line.strip().split('\t')
        cluster, contig = line[0],line[2]
        if contig not in contig_clusters:
            contig_clusters[contig] = cluster

### Parse through the fasta file to be subsetted and save the entries accordingly 
print('Parsing Contig Fasta File')
fasta_records = dict()
for record in SeqIO.parse(open(FASTA_FILE, 'r'), 'fasta'):
        cluster = None
        record.description = record.id
        contig = record.id
        if contig in contig_clusters:
            cluster = contig_clusters[contig]
            if cluster not in fasta_records:
                fasta_records[cluster] = []
                fasta_records[cluster] += [record]
            else:
                fasta_records[cluster] += [record]
contig_clusters = None

### Parse through fasta records of each cluster and write out to seperate bins
print('Writing fasta entries to seperate bins')
for cluster in fasta_records.keys():
    outhandle = open('{}/{}.fna'.format(OUTDIR,cluster),  'w')
    for record in fasta_records[cluster]:
        SeqIO.write(record, outhandle, 'fasta')