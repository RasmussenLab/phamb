#!/bin/python
import argparse
import numpy as np
import os
import sys
import pathlib
import csv
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-v', help='checkv directory')
parser.add_argument('-o', help='Viral proteome file')




def write_checkv_population_type(args):
    """[summary]

    Args:
        args ([type]): [description]
        bins ([type]): [description]
    """
    quality_file = os.path.join(args.v,'quality_summary.tsv')

    population_type_score = {3:'HQ-ref',2:'MQ-ref',1:'Dark-matter'}
    viral_population_types = dict()
    with open(quality_file,'r') as infile:
        reader = csv.DictReader(infile,delimiter='\t') 
        for row in reader:
            contigid  = row['contig_id']
            motherbin = contigid.split('_')[1]
            population_type = None

            ### If not a Cluster contains a single bin with at least MQ quality -> Dark-matter
            ### If a cluster contains at least 1 bin that's HQ by the AAI method -> HQ-ref 
            ### elif cluster contains at least 1 bin that's MQ by the AAI method -> MQ-ref 
            if not row['checkv_quality'] in ['Medium-quality','High-quality','Complete']:
                population_type = 1
                if motherbin in viral_population_types:
                    if population_type > viral_population_types[motherbin]:
                        viral_population_types[motherbin] = population_type
                else:
                    viral_population_types[motherbin] = population_type
            else:
                if row['completeness_method'] == 'AAI-based' and float(row['genome_copies']) < 1.25:
                    population_type = 2 
                    if row['checkv_quality'] in ['High-quality','Complete']:
                        population_type = 3
                    if motherbin in viral_population_types:
                        if population_type > viral_population_types[motherbin]:
                            viral_population_types[motherbin] = population_type
                    else:
                        viral_population_types[motherbin] = population_type
                
    fileout = os.path.join(args.v,'new_population_types.tsv')
    with open(fileout,'w') as out:
        for motherbin in viral_population_types:
            _type = viral_population_types[motherbin]
            population_type = population_type_score[_type]
            lineout = '\t'.join( [motherbin,population_type] )
            out.write(lineout + '\n')



def orf_distribution(sequence_regions,gene_coords,total_genes, viral_genes, host_genes):
    viral_orfs = 0
    host_orfs = 0

    if sequence_regions[0] == 'NA':
        host_orfs = host_genes
        viral_orfs = viral_genes
    else:
        for i, regiontype in enumerate(sequence_regions):
            if regiontype == 'viral':
                start, end = gene_coords[i].split('-')
                viral_orfs += int(end) - int(start)
            elif regiontype == 'host':
                start, end = gene_coords[i].split('-')
                host_orfs += int(end) - int(start)
    viral_percentage = round( (viral_orfs/total_genes)*100,2)
    host_percentage = round( (host_orfs/total_genes)*100, 2)
    NA_percentage = round( 100-viral_percentage-host_percentage,2)
    return (viral_percentage, host_percentage,NA_percentage)

def read_checkv_quality(args):
    '''
    1. Load annotation of each Bin from the checkV quality_summary file 
    2. Define Pro-viral Regions 
    '''

    quality_file = os.path.join(args.v,'quality_summary.tsv')
    #quality_file = '01_checkv/quality_summary.tsv'
    viruses = dict()
    with open(quality_file,'r') as infile:
        reader = csv.DictReader(infile,delimiter='\t') 
        for row in reader:
            contigid  = row['contig_id']
            if not row['checkv_quality'] in ['Medium-quality','High-quality','Complete']:
                continue
            viruses[contigid] = row
            
    ### Parse Contamination file and define gene coordinates for Viral Regions in each Contig
    ### There is a case 
    contamination_file = os.path.join(args.v,'contamination.tsv')
    with open(contamination_file,'r') as infile:
        reader = csv.DictReader(infile,delimiter='\t') 
        for row in reader:
            contigid  = row['contig_id']
            if not contigid in viruses:
                continue
            
            regions = row['region_types'].split(',')
            ### Formatting differences between most recent version of CheckV 
            if not 'region_genes' in row:
                gene_coords = row['region_coords_genes'].split(',')
                viruses[contigid]['viral_length'], viruses[contigid]['host_length']  = row['proviral_length'], row['host_length']
                circular = True if 'DTR' in viruses[contigid]['completeness_method'] else False
            else:
                gene_coords = row['region_genes'].split(',')
                viruses[contigid]['viral_length'], viruses[contigid]['host_length']  = row['viral_length'], row['host_length']
                total_genes = int(viruses[contigid]['gene_count'])
                circular = True if 'DTR' in viruses[contigid]['termini'] else False
            
            total_genes, viral_genes, host_genes = int(viruses[contigid]['gene_count']), int(row['viral_genes']), int(row['host_genes'])
            viruses[contigid]['viral_host_percentage'] = orf_distribution(regions, gene_coords,total_genes, viral_genes, host_genes )
            if regions[0] == 'NA':
                if int(viral_genes) > int(host_genes) or int(host_genes) == 0 or viruses[contigid]['completeness_method'] in ['AAI-based (high-confidence)','AAI-based']:
                    seq = np.arange(1, int(total_genes)+1) 
                    viruses[contigid]['viral_gene_coords'] = {v:k for k,v in enumerate(seq)}
            elif circular:      # Exception for DTR circular genomes 
                seq = np.arange(1, int(viruses[contigid]['gene_count'])+1) 
                viruses[contigid]['viral_gene_coords'] = {v:k for k,v in enumerate(seq)} 
            else:
                viral_regions_indices = []
                for i, regiontype in enumerate(regions):
                    if regiontype == 'viral':
                        viral_regions_indices.append(i)
                
                viral_gene_coords = np.array([])
                for i in viral_regions_indices:
                    start, end = gene_coords[i].split('-')
                    seq = np.arange(int(start),int(end)+1)
                    viral_gene_coords = np.append(viral_gene_coords,seq)
                viruses[contigid]['viral_gene_coords'] =  {v:k for k,v in enumerate(viral_gene_coords)}

    return viruses

def write_viral_proteomes(args,viruses):
    '''
    Parse Proteins predicted and write out proteins for MQ>= Viruses
    '''
    virus_proteins = os.path.join(args.o)
    proteins = os.path.join(args.v,'tmp','proteins.faa')
    ### Write out proteins 
    with open(virus_proteins, "w") as handle:
        for record in SeqIO.parse(proteins, "fasta"):
            header = record.id
            gene_number = header.split('_')[-1]
            contigid = '_'.join( header.split('_')[:-1] )
            if not contigid in viruses:
                continue
            if 'viral_gene_coords' in viruses[contigid]:
                if int(gene_number) in viruses[contigid]['viral_gene_coords']:
                    SeqIO.write(record, handle, "fasta")

def write_viral_contigs(args,viruses):
    '''
    Parse Cleaned contigs made by CheckV - subset to contigs from MQ >= Viruses
    '''
    checkv_cleaned_contigs = os.path.join(args.v,'cleaned_contigs.fna')
    virus_contigs = os.path.join(args.v,'virus_contigs.fna')
    ### Write out proteins 
    with open(virus_contigs, "w") as handle:
        for record in SeqIO.parse(checkv_cleaned_contigs, "fasta"):
            contigid = record.id 
            if not contigid in viruses:
                continue
            SeqIO.write(record, handle, "fasta")


def write_viral_host_contamination(args,viruses):

    fileout = os.path.join(args.v,'viral_host_contamination.tsv')
    with open(fileout,'w') as out:
        header = ['contig_id','contig_length','viral_length','host_length','perc_viral','perc_host','perc_NA','method','checkv_quality']
        out.write('\t'.join(header) + '\n')
        for contigid in viruses:
            virus = viruses[contigid]
            method = virus['completeness_method'].split(' ')[0]
            quality = virus['checkv_quality']
            contig_length, viral_length, host_length =  virus['contig_length'], virus['viral_length'], virus['host_length']
            viral_host_percentage = viruses[contigid]['viral_host_percentage']
            perc_viral, perc_host, perc_NA = viral_host_percentage[0], viral_host_percentage[1], viral_host_percentage[2]
            if viral_length == 'NA':
                viral_length, host_length = round(int(contig_length)*(perc_viral/100),2) ,round(int(contig_length)*(perc_host/100),2) 
            lineout = [contigid,contig_length, viral_length, host_length, perc_viral, perc_host, perc_NA,method,quality]
            lineout = [str(x) for x in lineout]
            out.write('\t'.join(lineout)+'\n')

                 

if __name__ == "__main__":
    '''
    1. Identify and organise Quality Viruses and Proviruses
    2. Write out their proteomes to a seperate file for annotation with I.e. Interproscan or Taxonomy.  
    '''
    args = parser.parse_args()
    write_checkv_population_type(args)

    viruses = read_checkv_quality(args)
    write_viral_proteomes(args,viruses)
    write_viral_contigs(args,viruses)
    write_viral_host_contamination(args,viruses)






class meh:
    def __init__(self,name):
        self.v = name 

args = meh('07_binannotation/checkv/VAMB_bins')




