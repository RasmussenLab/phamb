#!/bin/python
import argparse
import os
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-t',help='Table with Cluster identifiers etc.')
parser.add_argument('-i',help='CheckV results directory')
parser.add_argument('-o',help='Out Directory for .fna files')



def read_cluster_info(args):
    '''
    1. Loads a table with a Cluster identifier for all novo Clusters 
    2. Loads the genome entries into a dict based on the list of Clusters from the CheckV file. 
    '''

    mbins = set()
    with open(args.t,'r') as infile:
        header = infile.readline()

        for line in infile:
            motherbin, ngenomes, median_size , mad, mad_perc = line.strip().split('\t')
            if float(mad_perc)<20 and int(ngenomes)>2:
                mbins.add(motherbin)
    
    ### Parse CheckV standard file 
    quality_bins = {k:set() for k in mbins}
    checkv_quality = os.path.join(args.i,'quality_summary.tsv')
    with open(checkv_quality,'r') as infile:
        header = infile.readline().strip().split('\t')
        genome_copies_index = header.index('genome_copies')
        contamination_index = header.index('contamination')
        quality_index = header.index('checkv_quality')
        method_index = header.index('completeness_method')
        miuvig_quality = header.index('miuvig_quality')
        viral_genes_index = header.index('viral_genes')
        prophage_index = header.index('prophage')
        completeness_index = header.index('completeness')
        
        for line in infile:
            line = line.strip().split('\t')
            binid = line[0]
            
            motherbin = binid.split('_')[1]

            if not motherbin in mbins:
                continue
                
            quality = line[quality_index]
            completeness = line[completeness_index]
            genome_copies = float(line[genome_copies_index])
            contaminaton = line[contamination_index]

            ### Filter off Bins with little Viral evidence
            if quality == 'Not-determined' and line[viral_genes_index] == '0':
                continue
            if completeness =='NA':
                continue

            if float(completeness) > 75 and genome_copies < 1.25:
                quality_bins[motherbin].add(binid)
    return quality_bins

def write_cluster_fasta(args,quality_bins):
    '''
    Parse cleaned_contiggs.fna and write HQ-Bins to .fna file 
    
    '''


    ### Load Entries from cleaned_contigs.fna 
    ### Remember to choose a representative genome (largest one of the cluster)
    bins_fasta = os.path.join(args.i,'cleaned_contigs.fna')

    bins_dict = {}
    for record in SeqIO.parse(open(bins_fasta, 'r'), 'fasta'):
        recordname = record.id
        motherbin = recordname.split('_')[1]
        if motherbin in quality_bins:
                if recordname in quality_bins[motherbin]:
                    if not motherbin in bins_dict:
                        bins_dict[motherbin] = {}
                        sequence_len = len(record.seq)
                        bins_dict[motherbin]['repgenome'] = (recordname,sequence_len)
                        bins_dict[motherbin][recordname] = record
                    else:
                        sequence_len = len(record.seq)
                        sequence_len_rep = bins_dict[motherbin]['repgenome'][1]
                        if sequence_len > sequence_len_rep:
                            bins_dict[motherbin]['repgenome'] = (recordname,sequence_len)
                        bins_dict[motherbin][recordname] = record

    ### Write out Entries to seperate .fna files
    if not os.path.exists(args.o): 
        os.mkdir(args.o)
    for motherbin in bins_dict:
        motherbinfile = os.path.join(args.o,motherbin+'.fna')
        motherbinfile_rep = os.path.join(args.o,motherbin+'.rep.fna')

        with open(motherbinfile_rep , 'w') as outhandle:
            recordname = bins_dict[motherbin]['repgenome'][0]
            record = bins_dict[motherbin][recordname]
            SeqIO.write(record, outhandle, 'fasta')

        with open(motherbinfile , 'w') as outhandle:
            for recordname in bins_dict[motherbin]:
                if recordname != 'repgenome':
                    record = bins_dict[motherbin][recordname]
                    SeqIO.write(record, outhandle, 'fasta')




if __name__ == "__main__":
    args = parser.parse_args()

    quality_bins = read_cluster_info(args)
    write_cluster_fasta(args,quality_bins)


