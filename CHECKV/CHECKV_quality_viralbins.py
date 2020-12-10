#!/bin/python

import argparse
import os
import sys
import pathlib
import csv
from Bio import SeqIO
import gzip
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-c',help='VAMB directory')
parser.add_argument('-a',help='Annotation summary files directory')
parser.add_argument('-v', help='checkv directory')


class contig_class:
    def __init__(self):
        pass

class genome_class:
     def __init__(self,bin_name):
        self.bin = bin_name
        self.viral_length = 0
        self.contamination = "NA"
        self.total_genes = "NA"
        self.viral_genes = "NA"
        self.host_genes = "NA"
        self.completeness = "NA"
        self.complete = False
        self.method = "NA"
        

class Args:
    def __init__(self,v,c,o,db):
        self.v = v 
        self.o = o
        self.c = c
        self.db = db


args = Args('07_binannotation/checkv/VAMB_contigs','05_binning/vamb_on_jgi_v3/HMP2','lel','/home/projects/cpr_10006/projects/phamb/databases/checkv/checkv-db-v0.6')


if __name__ == "__main__":
    args = parser.parse_args()

    ### Load additional info about each Genome 
    contig_contamination = dict()
    p = os.path.join(args.v, "contamination.tsv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            contigid = r["contig_id"]
            contig = contig_class()
            contig.viral_length = int(r["viral_length"])
            contig.viral_genes = int(r["viral_genes"])
            contig_contamination[contigid] = contig
    
    p = os.path.join(args.v, "VAMB_completeness_all.tsv")
    VAMB_bins = dict()
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            genome = genome_class(r["vamb_bin"])
            contigs = r["contigs_in_genome"].split(';')
            viral_length, viral_genes = 0, 0
            for c in contigs:
                viral_genes += contig_contamination[c].viral_genes
                viral_length += contig_contamination[c].viral_length
        
            genome.viral_length = viral_length
            genome.viral_genes = viral_genes

            ### Continue if already circular entry noted.
            if genome.bin in VAMB_bins:
                continue 

            # AAI-based esimtate (med/high-confidence)
            if r["aai_confidence"] in ("high", "medium"):
                genome.completeness = round(float(r["aai_completeness"]), 2)
                genome.method = "AAI-based (%s-confidence)" % r["aai_confidence"]
            # HMM-based estimate
            elif r["hmm_lower_completeness"] != "NA":
                genome.completeness = round(float(r["hmm_lower_completeness"]), 2)
                genome.method = "HMM-based (lower-bound)"
            # no completeness estimate whatsoever
            else:
                genome.completeness = "NA"
                genome.method = "NA"

            # Is the genome circular?
            if '_circular' in r["vamb_split"]:
                genome.method = 'circular_genome'
                genome.completeness = 100
                genome.complete = True
            
            VAMB_bins[genome.bin] = genome

    for VAMB_bin in VAMB_bins:
        genome = VAMB_bins[VAMB_bin]
        if genome.complete:
            genome.quality = "Complete"
            genome.miuvig = "High-quality"
        elif genome.completeness == "NA":
            genome.quality = "Not-determined"
            genome.miuvig = "Genome-fragment"
        elif genome.completeness >= 90:
            genome.quality = "High-quality"
            genome.miuvig = "High-quality"
        elif genome.completeness >= 50:
            genome.quality = "Medium-quality"
            genome.miuvig = "Genome-fragment"
        elif genome.completeness < 50:
            genome.quality = "Low-quality"
            genome.miuvig = "Genome-fragment"
    
    header = [
        "vamb_bin",
        "viral_length",
        "viral_genes",
        "checkv_quality",
        "miuvig_quality",
        "completeness",
        "completeness_method"
    ]
    outfile = os.path.join(args.v,'VAMB_quality_summary_all.tsv')
    with open(outfile,'w') as out:
        out.write("\t".join(header) + "\n")
        for VAMB_bin in VAMB_bins:
            genome = VAMB_bins[VAMB_bin]
            row = [
                genome.bin,
                genome.viral_length,
                genome.viral_genes,
                genome.quality,
                genome.miuvig,
                genome.completeness,
                genome.method
            ]
            out.write("\t".join([str(_) for _ in row]) + "\n")
        