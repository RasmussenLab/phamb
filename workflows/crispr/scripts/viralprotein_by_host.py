#!/bin/bash

import argparse
import os
import sys
import csv
from collections import defaultdict
import numpy as np

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-v',help='CheckV directory')
parser.add_argument('-g', help='GTDB-TK taxonomy file')
parser.add_argument('-p',help='Prophage/Host annotation directory')



def load_MAG_taxonomy(args):
    '''
    GTDB-TK annotation of NC MAGs.
    IDEALLY this should be a consensus vote, not the first occuring...
    If it's annotated all the way to Species level don't parse anymore.
    '''

    gtdb_file = os.path.join(args.g,'mag_taxonomy.txt')
    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    #dlineage = {lin:None for lin in lineage}

    parsed_motherbins = set()
    MAG_tax = dict()
    with open(gtdb_file,'r') as infile:
        header = infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            motherbin = line[0].split('_')[1]

            if motherbin in MAG_tax:
                if MAG_tax[motherbin]['species'] != 'NA':
                    continue

            MAG_tax[motherbin] = {lin:None for lin in lineage}
            tax = line[2]
            tax = tax.split(';')
            tax = [ x.split('__')[1] for x in tax ]
            for i,lin in enumerate(lineage):
                entry = tax[i]
                if entry == '':
                    entry = 'NA'
                if '_' in entry:
                    entry = entry.split('_')[0]
                if entry == 'Bacteroidota':
                    entry = 'Bacteroidetes'
                
                MAG_tax[motherbin][lin] = entry
    return MAG_tax




class viral_proteome:
     def __init__(self,motherbin):
        self.motherbin = motherbin
        self.MAGs = set()
        self.viral_proteins = dict()


def load_viral_host_annotation(args):
    '''
    1. Load MAG-taxonomy 
    2. Connect Viruses to MAGs 
    '''
    MAG_tax = load_MAG_taxonomy(args)
    #mag_genera = dict() 
    viral_hosts = dict()
    filein = os.path.join(args.p , 'crispr_prophage_results' ,'MAG_crispr_prophage_summary.5000.txt')
    
    ### Connect Viruses to MAGs
    with open(filein,'r') as infile:
        for line in infile:
            MAG_motherbin, MAG, viral_motherbin, viral_type, host_annotation, coverage = line.strip().split('\t')
            if viral_type in ['HQ-ref'] and host_annotation in ['both','prophage']:
                if not MAG_motherbin in MAG_tax:
                    continue
                #genus = MAG_tax[MAG_tax]['genus']
                
                if not viral_motherbin in viral_hosts:
                    VP = viral_proteome(viral_motherbin)
                    VP.MAGs.add(MAG_motherbin)
                    viral_hosts[viral_motherbin] = VP
                else:
                    viral_hosts[viral_motherbin].MAGs.add(MAG_motherbin)

    return viral_hosts
        
def parse_viral_protein_annotations(args, viral_hosts):
    '''

    '''
    MAG_tax = load_MAG_taxonomy(args)

    files = [os.path.join(args.v,'proteome_annotation',x) for x in ['MQLQND.bins.interpros.txt','nc.bins.interpros.txt']]
    for f in files:
        for r in csv.DictReader(open(f), delimiter="\t"):
            dbentry = r['db']
            dbdesc = r['dbid']
            viral_motherbin = r['motherbin']

            if not viral_motherbin in viral_hosts:
                continue
            VP = viral_hosts[viral_motherbin]
            VP.viral_proteins[dbentry] = dbdesc
    return viral_hosts


class host_taxonomy:
    def __init__(self,lineage):
        self.lineage = {lin:lineage[lin] for lin in ['superkingdom','phylum','class','order','family','genus','species'] }
        self.MAGs = set()
        self.pfam_count = defaultdict(int)
            
def count_by_taxonomic_host_rank(args, viral_hosts,taxonomic_rank = 'genus'):
    '''

    '''

    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    ### Load MAG taxonomy 
    MAG_tax = load_MAG_taxonomy(args)
    all_pfams = set()
    host_dict = dict()
    for MAG_motherbin in MAG_tax:
        MAG_lineage =  MAG_tax[MAG_motherbin]
        genera = MAG_lineage[taxonomic_rank]

        if not genera in host_dict:
            H = host_taxonomy(MAG_lineage)
            H.MAGs.add(MAG_motherbin)
            host_dict[genera] = H
        else:
            H = host_dict[genera]
            H.MAGs.add(MAG_motherbin)
        
    ### Which Taxa does the Virus infect  
    for viral_motherbin,VP in viral_hosts.items():
        pfams = VP.viral_proteins
        MAGs = VP.MAGs
        MAG_genera = set([MAG_tax[M][taxonomic_rank] for M in MAGs])

        ### Increment the PFAM count for each Genera associated with the Virus  
        for gen in MAG_genera:
            for pf in pfams:
                all_pfams.add(pf)
                if gen in host_dict:
                    host_dict[gen].pfam_count[pf] += 1 
    
    return host_dict, all_pfams
    

def write_count_matrix(args,all_pfams,host_dict):


    ### 
    nentries = len(all_pfams)
    ntaxa = len(host_dict)
    taxa = list(host_dict)
    mat = np.zeros((nentries,ntaxa))

    long_table_list = []
    for i,pfam in enumerate(list(all_pfams)):
        for j,tax in enumerate(host_dict):
            mat[i,j] = host_dict[tax].pfam_count[pfam]
            long_table_list.append([pfam,tax, host_dict[tax].pfam_count[pfam]])



    rows = np.array( list(all_pfams) ) 
    rows = rows.reshape(-1,1)
    taxa.insert(0,'feature')
    header = np.array(taxa).reshape(1,-1)
    mat = np.hstack( (rows, mat))
    mat = np.vstack( (header, mat))
    fileout = 'test.matrix.txt'
    np.savetxt(fileout, mat, delimiter='\t', header= '',  newline='\n', fmt='%s')


        




class Args:
    def __init__(self,v,g,p):
        self.v = v 
        self.g = g
        self.p = p



if __name__ == "__main__":
    '''
    1. Establish connections between MAGs and viruses.
    2. Parse Viral Protein annotations (from HQ-viruses) associated with each MAG 
    3. Count up the number of different PFAMs by MAG-Genus - making a PFAM x MAG table. 
        - Count only a PFAM once pr. Viral Motherbin!
    '''
    
    
    args = parser.parse_args()

    args = Args('07_binannotation/checkv/VAMB_bins','07_binannotation/bacteria/gtdbtk','08_crisprcas')

    ### Establish connections between MAGs and viruses
    viral_hosts = load_viral_host_annotation(args)

    ### Parse Viral Protein annotations (from HQ-viruses) associated with each MAG 
    viral_hosts = parse_viral_protein_annotations(args ,viral_hosts)

    ### Count up the number of Distinct PFAMS bby MAG-tax
    host_dict,all_pfams  = count_by_taxonomic_host_rank(args, viral_hosts,taxonomic_rank = 'genus')
