#!/bin/python
import argparse
import os
import sys
from Bio import SeqIO
import subprocess
import numpy as np


parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-i',help='CheckV results directory')
parser.add_argument('-o',help='Out Directory for .fna files')

def write_all_viral_bins(args):
    '''
    Parse cleaned_contiggs.fna and write bins to .fna
    
    '''
    bins_fasta = os.path.join(args.i,'cleaned_contigs.fna')
    quality_file = os.path.join(args.i,'quality_summary.tsv')

    ### Determine the longest representative viral genome of each Cluster
    bins_representative = dict()
    with open(quality_file,'r') as infile:
         header = infile.readline()
         for line in infile:
            line = line.strip().split('\t')
            binid = line[0]
            motherbin = binid.split('_')[1]
            if not motherbin in bins_representative:
                bins_representative[motherbin] = set([binid])
            else:
                bins_representative[motherbin].add(binid)

    ### Parse fasta with cleaned contigs
    bins_dict = {}
    for record in SeqIO.parse(open(bins_fasta, 'r'), 'fasta'):
        recordname = record.id
        motherbin = recordname.split('_')[1]
        bin = recordname
        if motherbin in bins_representative:
                if bin in bins_representative[motherbin]:
                        if not motherbin in bins_dict:
                            bins_dict[motherbin] = {}
                            bins_dict[motherbin][bin] = record
                        else:
                            bins_dict[motherbin][bin] = record

    ### Write out Entries to seperate .fna files
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    directory_out = args.o
    viral_bins_file = os.path.join(directory_out,'all.viral.genomes.txt')
    with open(viral_bins_file,'w') as outfile:
            for motherbin in bins_dict:
                    clusterdirectory = os.path.join(args.o,motherbin)
                    if not os.path.exists(clusterdirectory):
                            os.mkdir(clusterdirectory)
                    for bin in bins_dict[motherbin]:
                            binfile = os.path.join(clusterdirectory,bin+'.fna')
                            if not os.path.exists(binfile):
                                    with open(binfile, 'w') as outhandle:
                                        record = bins_dict[motherbin][bin]
                                        SeqIO.write(record, outhandle, 'fasta')
                            outfile.write(binfile+'\n')


if __name__ == "__main__":
    args = parser.parse_args()
    write_all_viral_bins(args)