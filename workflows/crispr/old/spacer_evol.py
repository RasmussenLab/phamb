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

parser.add_argument('-m',help='Metadata file')
parser.add_argument('-t',help='Directory with CRISPR cluster matrices.')
parser.add_argument('-i',help='Directory with CRISPR spacers')
parser.add_argument('-o',help='Out Directory for .fna files')


def read_metadata(args):

    metadata = {}
    with open(args.m,'r') as infile:
        header = infile.readline()
        for line in infile:
            project, Sample, Subject, week_num, diagnosis, datatype = line.strip().split('\t')
            metadata[Sample] = [Subject,week_num,diagnosis]
    return metadata


def read_clustermatrices(args,bacterial_clusters):
    
    bacterial_prophage = dict() 
    for bacterial_motherbin in bacterial_clusters:
        bacterial_prophage[bacterial_motherbin] = dict()

    for bacterial_motherbin in bacterial_clusters:
        matrixfile = os.path.join(args.t,bacterial_motherbin+'.clustermatrix.prophage.txt')
        
        prophage_state = {0:'Nohit',1:'crisprhit',0.5:'prophage',1.5:'crisprhit_prophage'}

        header = None
        with open(matrixfile,'r') as infile:
            header = infile.readline().strip()
            bacterialbins = header.split('\t')[1:]
            samples = [b.split('_')[0] for b in bacterialbins]
            for sample in samples:
                bacterial_prophage[bacterial_motherbin][sample] = {}

            for line in infile:
                viral_cluster = line.strip().split('\t')[0]
                hits = line.strip().split('\t')[1:]
                hits = np.array(hits,dtype='float64')

                for i,sample in enumerate(samples):
                    state = prophage_state[hits[i]]
                    bacterial_prophage[bacterial_motherbin][sample][viral_cluster] = state
    return bacterial_prophage
        

def read_spacer_dict(args,bacterial_clusters):
    """[Reads a .fna file with spacers and maps the spacer-sequence to the entry naame ]

    Args:
        args ([type]): [description]
    """
    magspacerdict = dict()
    spacers_in_MAG = dict()
    spacerclstr = dict()

    for MAG in bacterial_clusters:
        mag_spacer_file = os.path.join(args.i,'spacer_by_MAG',MAG,MAG+'.spacers.fna')
        spacers_in_MAG[MAG] = dict()
        spacerid = None
        with open(mag_spacer_file,'r') as infile:
            for line in infile:
                if line[0] == '>':
                    spacerid = line[1:].strip()
                    bin = spacerid.split(':')[0]
                    sample =  bin.split('_')[0]
                    motherbin = bin.split('_')[1]
                    if not sample in spacers_in_MAG[motherbin]:
                        spacers_in_MAG[motherbin][sample] = [spacerid]
                    else:
                        spacers_in_MAG[motherbin][sample] += [spacerid]
                    
                else:
                    seq = line.strip()
                    magspacerdict[spacerid] = seq
        mag_spacer_clstr_file = os.path.join(args.i,'spacer_by_MAG',MAG,MAG+'.cdhit.clstr')

        spacerclstr[MAG] = dict()
        i = 0
        clstr = None
        with open(mag_spacer_clstr_file,'r') as infile:
            for line in infile:
                if line[0] == '>':
                    i += 1 
                    clstr = str(i)
                else:
                    line = line.strip().split('>')[1]
                    spacerid = line.split('...')[0]
                    spacerclstr[MAG][spacerid] = clstr

    return magspacerdict, spacers_in_MAG, spacerclstr



if __name__ == "__main__":
    args = parser.parse_args()

    bacterial_clusters = ['146']
    metadata = read_metadata(args)
    magspacerdict, spacers_in_MAG, spacerclstr = read_spacer_dict(args,bacterial_clusters)
    bacterial_prophage = read_clustermatrices(args,bacterial_clusters)

    outfile = 'subject_mag_spacsers.tsv'
    with open(outfile,'w') as out:

        for bacterial_motherbin in bacterial_clusters:
            for Sample in bacterial_prophage[bacterial_motherbin]:
                subject = metadata[Sample][0]
                week = metadata[Sample][1]

                for viral_motherbin in bacterial_prophage[bacterial_motherbin][Sample]:
                    state = bacterial_prophage[bacterial_motherbin][Sample][viral_motherbin]

                    spacer_clstrs = set()
                    for spacerid in spacers_in_MAG[bacterial_motherbin][Sample]:
                        if state != 'Nohit':
                            seq = magspacerdict[spacerid]
                            clstr = spacerclstr[bacterial_motherbin][spacerid]
                            spacer_clstrs.add(int(clstr))  
                            lineout = [subject,Sample,week,bacterial_motherbin,viral_motherbin,state,clstr]
                            out.write('\t'.join(lineout) +'\n')









