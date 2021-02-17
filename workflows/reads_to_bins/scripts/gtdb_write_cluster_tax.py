#!/bin/python

import sys
import os
import argparse

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
   ''')
parser.add_argument('-g', help='GTDB-TK taxonomy file')
parser.add_argument('-o', help='Taxonomy file for each VAMB cluster')


def load_MAG_taxonomy(args):
    '''
    GTDB-TK annotation of NC MAGs.
    IDEALLY this should be a consensus vote, not the first occuring...
    If it's annotated all the way to Species level don't parse anymore.
    '''

    gtdb_file = args.g
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

def write_taxonomy_file(args,MAG_tax):
    '''
    Write out taxonomy from GTDBTK for each Mothercluster
    '''
    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    fileout = args.o 

    with open(fileout,'w') as out:
        header = ['motherbin'] 
        header += lineage
        header = '\t'.join(header)+'\n'
        out.write(header)
        for MAG in MAG_tax:
            MAG_lineage = []
            for lin in lineage:
                MAG_lineage += [MAG_tax[MAG][lin]] 
            lineout = [str(MAG)] 
            lineout += MAG_lineage
            out.write('\t'.join(lineout)+'\n')

if __name__ == "__main__":
    args = parser.parse_args()
    MAG_tax = load_MAG_taxonomy(args)
    write_taxonomy_file(args,MAG_tax)



