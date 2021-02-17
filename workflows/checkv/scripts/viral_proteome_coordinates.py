#!/bin/python
import argparse
import os
import sys
import pathlib


parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
   ''')
parser.add_argument('-v', help='CheckV directory')
parser.add_argument('-d', help='CheckV Db directory')



class genome_annotation():
    def __init__(self,name):
        self.binname = name
        self.ORFs = []
        self.ORF_coordinates = dict()
        self.ORF_annotation = dict()

def load_checkvhmm_annotation(args):

    checkv_hmms = os.path.join(args.d,'hmm_db','checkv_hmms.tsv')
    checkv_hmms = '/home/projects/cpr_10006/projects/phamb/databases/checkv/checkv-db-v0.6/hmm_db/checkv_hmms.tsv'

    hmmlookup = dict()
    with open(checkv_hmms,'r') as infile:
        header = infile.readline().strip().split('\t')
        acc_index, hmm_index, desc_index  = header.index('acc'),header.index('hmm'), header.index('desc')
        for line in infile:
            line = line.strip().split('\t')
            acc, hmm, description = line[acc_index],line[hmm_index],line[desc_index]
            hmmlookup[hmm] = (acc,description)
    return hmmlookup



def write_cluster_coordinates(viral_clusters):
    '''
    1. Parse gene features (coordinates) and annotation (hmm hit)
    2. Store in Class object pr. Genome
    ''' 
    genomes = dict()
    gene_feature_file = os.path.join('07_binannotation/checkv/VAMB_bins','tmp','gene_features.tsv')
    gene_annotation_file = os.path.join('07_binannotation/checkv/VAMB_bins','tmp','gene_annotations.tsv')

    ### Parse Gene feeatures - coordinates and strand orientation
    with open(gene_feature_file,'r') as infile:
        header = infile.readline().strip().split('\t')
        bin_index, gene_index = header.index('contig_id'),header.index('gene_num')
        gene_start, gene_end, gene_strand = header.index('start'), header.index('end'), header.index('strand')
        for line in infile:
            line = line.strip().split('\t')
            binname, gene, start, end, strand  = line[bin_index], line[gene_index], line[gene_start], line[gene_end], line[gene_strand]

            cluster = binname.split('_')[1]
            if not cluster in viral_clusters:
                continue

            if not binname in genomes:
                genome = genome_annotation(binname)
            else:
                genome = genomes[binname]
            genome.ORFs += [gene]
            genome.ORF_coordinates[gene] = (start, end, strand)
            genomes[binname] = genome

    ### Parse Gene Annotations 
    with open(gene_annotation_file,'r') as infile:
        header = infile.readline().strip().split('\t')
        bin_index, gene_index = header.index('contig_id'),header.index('gene_num')
        hmm_cat, hmm_target = header.index('hmm_cat'),header.index('target_hmm')
        
        for line in infile:
            line = line.strip().split('\t')
            binname, gene, cat, target = line[bin_index], line[gene_index], line[hmm_cat], line[hmm_target]

            cluster = binname.split('_')[1]
            if not cluster in viral_clusters:
                continue
            if binname in genomes:
                genome = genomes[binname]
                genome.ORF_annotation[gene] = (cat,target)
                genomes[binname] = genome
    return genomes

def write_VMAG_proteomes(hmmlookup, viral_clusters,genomes):
    '''
    Write out a Table pr. Viral cluster with coordinates of each Protein and annotation
    '''

    ### 
    outdirectory = os.path.join('07_binannotation/checkv/VAMB_bins','viral_protein_coords')
    if not os.path.exists(outdirectory):
        os.makedirs(outdirectory)
        
    for cluster in viral_clusters:
        outfile = os.path.join('07_binannotation/checkv/VAMB_bins','viral_protein_coords',cluster + '_coords.txt')
        with open(outfile,'w') as out:
            header = '\t'.join(['genome','genenumber','start','end','orientation','gene','checkvcat','acc','description'])
            out.write(header+'\n')
            for binname in genomes:
                clustern = binname.split('_')[1]
                if clustern == cluster:
                    genome = genomes[binname]

                    for gene in genome.ORFs:
                        start, end, strand = genome.ORF_coordinates[gene]
                        cat, target, acc, desc = 'NA', 'NA', 'NA', 'NA'
                        if gene in genome.ORF_annotation:
                            cat, target = genome.ORF_annotation[gene]
                        if target in hmmlookup:
                            acc, desc = hmmlookup[target]
                        lineout = '\t'.join([binname,gene, start, end, strand, target, cat, acc, desc])
                        out.write(lineout+'\n')


if __name__ == "__main__":
    args = parser.parse_args()
    viral_clusters = ['653','502']

    hmmlookup = load_checkvhmm_annotation(args)
    genomes = write_cluster_coordinates(viral_clusters)
    write_VMAG_proteomes(hmmlookup,viral_clusters,genomes)





