#!/bin/python
import argparse
import os
import sys
import numpy as np
import collections
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
   ''')
parser.add_argument('-m', help='metadata file',default=None)


'''
Create Table with MAG subtypes 
'''

def parse_metadata():
    
    sample_meta = dict()

    with open('metadata','r') as infile:
        header = infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            sample, subject, week, diagnosis = line[1],line[2],line[3],line[4]
            sample_meta[sample] = '_'.join([subject, week, diagnosis])
    return sample_meta
            



def get_MAG_overview():
    """[Read CheckM file of MQNC MAGs]

    Args:
        args ([type]): [description]
    """
    MAGs = dict()
    checkm_file = '07_binannotation/bacteria/MQNC_MAGS.txt'
    with open(checkm_file,'r') as infile:
        for line in infile:
            line = line.strip()
            genomeid = os.path.basename(line)
            genomeid = genomeid.replace('.fna','')
            MAGs[genomeid] = dict()

    return MAGs

class cluster_crispr():
    def __init__(self,name):
        self.clustername = name
        self.subtypes = []
        self.binsubtypes = dict()
        self.spacers = dict()
        self.majority_vote = None


def parse_subtypes(args,MAGs):
    '''
    1. Read predicted CAS subtypes for all MAGs 
    2. Compute subtype frequencies and Subtype Majority vote for each MAG cluster
    '''

    ### load metadat 
    if not args is None: 
        sample_meta = parse_metadata()
    else:
        sample_meta = None

    ### Load spacers
    spacers = dict()
    spacerfasta = os.path.join('07_binannotation/bacteria/snakecrisprcasfinder','all.spacers.fna')
    for record in SeqIO.parse(spacerfasta, "fasta"):
        header = record.id 
        contig = header.split(':')[0]
        if not contig in spacers:
            spacers[contig] = dict()
            spacers[contig][header] = str(record.seq)
        else:
            spacers[contig][header] = str(record.seq)

    cluster_subtypes = dict()
    linesout = []
    ### Parse predictions
    for MAG in MAGs:
        sample, vamb_cluster = MAG.split('_')
        if not vamb_cluster == '216':
            continue        
        sample_metadata  = 'None'
        if not sample_meta is None:
            if sample in sample_meta:
                sample_metadata = sample_meta[sample]
        crispr_cas_file = os.path.join('07_binannotation/bacteria/snakecrisprcasfinder',MAG,'CRISPR_Cas.tab')
        if os.path.exists(crispr_cas_file):
            with open(crispr_cas_file,'r') as infile:
                header = infile.readline()
                for line in infile:
                    line = line.strip().split('\t')
                    contig, cas_prediction, crisprfile = line[0] , line[3], line[4]
                    crisprs = crisprfile.replace('[','').replace(']','').replace("'",'').replace(' ','').split(',')

                    subject, week, diagnosis = sample_metadata.split('_')
                    if not vamb_cluster in cluster_subtypes:
                        cc = cluster_crispr(vamb_cluster)
                        cc.subtypes += [cas_prediction]
                        cc.binsubtypes[MAG] = (cas_prediction, sample_metadata )
                        cluster_subtypes[vamb_cluster] = cc
                    else:
                        cluster_subtypes[vamb_cluster].subtypes += [cas_prediction]
                        cluster_subtypes[vamb_cluster].binsubtypes[MAG] = (cas_prediction, sample_metadata)

                    ### Get spacers 
                    for c in crisprs:
                        for s in spacers[c]:
                            cluster_subtypes[vamb_cluster].spacers[s] = spacers[c][s]
                            line = [MAG,subject,week,diagnosis,cas_prediction,s, spacers[c][s] ]
                            linesout += [line]
    
    spacer_file_out = os.path.join('07_binannotation/bacteria/snakecrisprcasfinder','216.bvulgatus.spacers.txt')
    with open(spacer_file_out,'w') as out:
        out.write('\t'.join(['MAG','subject','week','diagnosis','cassubtype','contigid','spacer'])+'\n')
        for l in linesout:
            out.write('\t'.join(l)+'\n')
    
    ### Calculate frequencies and majority vote for each Cluster
    for vamb_cluster in cluster_subtypes:
        cc =  cluster_subtypes[vamb_cluster]
        subtypes = cc.subtypes
        subtypes_frequencies = collections.Counter(subtypes)
        subtypes_frequencies = {k: v/len(subtypes) for k,v in subtypes_frequencies.items()}
        sorted_subtypes = [k for k, v in sorted(subtypes_frequencies.items(), key=lambda item: item[1], reverse=True)]
        #print(vamb_cluster,'...',subtypes_frequencies, '...' )
        cluster_subtypes[vamb_cluster].majority_vote = sorted_subtypes[0]

    return cluster_subtypes
    

def write_subtypes(cluster_subtypes):

    fileout = os.path.join('08_crisprcas','crispr_prophage_results','cctyper_subtypes.txt')

    with open(fileout,'w') as out:
        for vamb_cluster in cluster_subtypes:
            cc = cluster_subtypes[vamb_cluster]
            majority_type = cc.majority_vote
            binsubtypes = cc.binsubtypes
            for MAG in binsubtypes:
                cas_prediction, sample_metadata = binsubtypes[MAG]

                lineout = '\t'.join( [vamb_cluster, majority_type, MAG,cas_prediction, sample_metadata] )
                out.write(lineout+'\n')



if __name__ == "__main__":

    args = parser.parse_args()
    
    ### Create MAG contig overview 
    MAGs = get_MAG_overview()

    ### Parse CRISPR_cas.tab files of each MAG
    cluster_subtypes = parse_subtypes(args,MAGs)
    write_subtypes(cluster_subtypes)
