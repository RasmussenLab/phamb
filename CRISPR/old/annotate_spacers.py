#!/bin/python
import subprocess
import os
import sys
import argparse


parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
   ''')
parser.add_argument('-v', help='CheckV output directory')
parser.add_argument('-i', help='Directory with spacer-blast files')
parser.add_argument('-p', help='Prophage annotation file')


def read_spacer_dict(args):
    """[Reads a .fna file with spacers and maps the spacer-sequence to the entry naame ]

    Args:
        args ([type]): [description]
    """
    spacerfile = os.path.join(args.i, 'new.all.spacers.fna')
    spacers_in_MAG = dict()
    magspacerdict = dict()
    spacerid = None
    with open(spacerfile,'r') as infile:
        for line in infile:
            if line[0] == '>':
                spacerid = line[1:].strip()
                bin = spacerid.split(':')[0]
                motherbin = bin.split('_')[1]
                if not motherbin in spacers_in_MAG:
                    spacers_in_MAG[motherbin] = [spacerid]
                else:
                    spacers_in_MAG[motherbin] += [spacerid]
                
            else:
                seq = line.strip()
                magspacerdict[spacerid] = seq

    return magspacerdict



def read_cluster_viruses(args):
    """[summary]

    Args:
        args ([type]): [description]
    """

    ### Parse CheckV standard file 
    checkv_quality = os.path.join(args.v,'quality_summary.tsv')
    bins_qualities = dict()

    ### Parse Prophage annotations (A)
    prophage_annotations = args.p
    with open(prophage_annotations,'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            viral_motherbin = line[3]
            MAG = line[0]
            bins_qualities[viral_motherbin] = 'Prophage'
    
    ### Score = [1,2,3,4] = [HQ,MQ,LQ,Not-Determined]
    score_tiers = {1:'HQ',2:'MQ',3:'LQ',4:'Not-Determined',0:'Prophage'}
    quality_score_dict = dict()
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
            if motherbin in bins_qualities:
                continue
                
            quality = line[quality_index]
            completeness = line[completeness_index]
            genome_copies = float(line[genome_copies_index])
            contaminaton = line[contamination_index]
            
            score = None

            ### Filter off Bins with little Viral evidence
            if quality == 'Not-determined' or completeness == '0.0':
                score = 4
                quality_score_dict[motherbin] = score 
                continue

            if float(completeness) > 0:
                score = 3
            if float(completeness) >= 50:
                score = 2
            if float(completeness) >= 90:
                score = 1
            if not motherbin in quality_score_dict:
                quality_score_dict[motherbin] = score
            else:
                if score < quality_score_dict[motherbin]:
                    quality_score_dict[motherbin] = score 
    
    ### catalogue each motherbin 
    for motherbin in quality_score_dict:
        score = quality_score_dict[motherbin]
        tier = score_tiers[score]
        bins_qualities[motherbin] = tier
    
    return bins_qualities

def parse_CRISPR_hits(args):
    """Function for parsing CRISPR-spacer blast-files

    Args:
        args (arguments from commandline): [arguments from commandline]

    """
    ### Parse MAG spacers
    magspacerdict = read_spacer_dict(args)
    bins_qualities = read_cluster_viruses(args)

    ### Parse blastn file of SPACEROME
    
    ### We should also Organise In a dictionary what the Spacer-hit corresponds to? i.e.
    # Known Viral Genome 
    # A known viral prophage?
    # No signal in CheckV at all?
    # A known viral prophage?
    # Then these annotations can be split pr. MAG/Higher taxonomy afterwards.
    # Also do we have Spacers from the same MAG targetting a bin across time? 

    spacer_out_file = os.path.join(args.i,'new.spacers.annotated.txt')
    spaceromefile = os.path.join(args.i, 'new.spacers.viralfragment.m6')
    with open(spacer_out_file,'w') as outfile:
        for f in [spaceromefile]:
            with open(f,'r') as infile:
                for line in infile:
                    line = line.strip().split()
                    spacername = line[0]
                    MAG = spacername.split(':')[0]
                    MAGmotherbin = MAG.split('_')[1]
                    phagebin = line[1]
                    phagemotherbin = phagebin.split('_')[1]
                    identity = float(line[2])
                    spacer_covered,spacer_len = int(line[3]),int(line[-2])
                    spacer_coverage = round((spacer_covered/spacer_len)*100,2)
                    bitscore = float(line[-3])
                    evalue = float(line[-4])
                    mismatches = int(line[4])

                    evidence = spacername.split('_')[-1]
                    evidence = 'evidence_' + evidence

                    ### Stringent filtering of hits
                    quality_tier = bins_qualities[phagemotherbin]
                    seq = magspacerdict[spacername]
                    lineout = [MAGmotherbin,MAG,phagebin,spacername,evidence, quality_tier, str(identity),str(spacer_coverage),str(bitscore), seq]
                    outfile.write('\t'.join(lineout)+'\n')
                    

if __name__ == "__main__":
    args = parser.parse_args()

    ###
    parse_CRISPR_hits(args)