#!/bin/bash 
import vamb
import argparse

### Usage 
parser = argparse.ArgumentParser(description='''
   Takes a .clusters.tsv file from VAMB and subsets the fasta file with contigs used in VAMB
   ''')
parser.add_argument('-c', help='Cluster file')
parser.add_argument('-f', help='Combined assembly')
parser.add_argument('-d', help='Out Directory of fasta files', type=str, default=None)
parser.add_argument('-m', help='Minimum MAG size',default=1*10**6)




if __name__ == "__main__":
    args = parser.parse_args()
    fastafile = args.f
    minimum_genome_size = args.m
    ### Load Fastadict 
    with vamb.vambtools.Reader(fastafile, 'rb') as filehandle:
        fastadict = vamb.vambtools.loadfasta(filehandle)

    clusterfile = args.c
    with open(clusterfile,'r') as filehandle:
        bins = vamb.vambtools.read_clusters(filehandle ,min_size=1)

    ### Split clusters to seperate genomes
    bins_splitted = vamb.vambtools.binsplit(bins, '_')

    binsdirectory = args.d

    ### The maxbins is a strange argument, set it to unfathomably high. 
    vamb.vambtools.write_bins(binsdirectory, bins_splitted, fastadict, compressed=False, maxbins=1*10**7, minsize=1*10**6)
