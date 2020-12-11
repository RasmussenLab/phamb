#!/bin/python
import os
import sys

def read_write_MAG_taxonomy():
    '''
    GTDB-TK annotation of NC MAGs.
    '''

    gtdb_file = '07_binannotation/bacteria/nc_gtdbtk_MQNC/taxonomy.txt'
    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    #dlineage = {lin:None for lin in lineage}

    if not os.path.exists:
        print('No GTDBTK taxonomy available? Generate this file', gtdb_file)
        sys.exit(1)
    parsed_motherbins = set()
    MAG_tax = dict()
    with open(gtdb_file,'r') as infile:
        header = infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            genomeid = line[0]
            motherbin = genomeid.split('_')[1]
            if motherbin in parsed_motherbins:
                continue
            parsed_motherbins.add(motherbin)
            MAG_tax[motherbin] = {lin:None for lin in lineage}
            tax = line[2]
            
            tax = tax.split(';')
            tax = [ x.split('__')[1] for x in tax ]
            for i,lin in enumerate(lineage):
                entry = tax[i]
                if entry == '':
                    entry = 'NA'
                if entry == 'Bacteroidota':
                    entry = 'Bacteroidetes'
                
                MAG_tax[motherbin][lin] = entry
    fileout1 = '10_abundances/VAMBV3.MAGmotherbin.tax.txt'

    with open(fileout1,'w') as out:
        out.write('motherbin' + '\t' + '\t'.join(lineage) + '\n')
        for MAG in MAG_tax:
            taxstring = [ MAG_tax[MAG][lin] for lin in lineage ]
            out.write(MAG + '\t' + '\t'.join(taxstring) + '\n')

if __name__ == "__main__":
    read_write_MAG_taxonomy()
