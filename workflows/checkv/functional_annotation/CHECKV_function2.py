#!/bin/python
import argparse
import os
import sys
import pathlib
from Bio import SeqIO
from collections import defaultdict
import re
import numpy as np


parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')


parser.add_argument('-v', help='checkv directory')
parser.add_argument('-o', help='Directory for files out')



def parse_protein_clstrs(args):
    """[summary]

    Args:
        args ([type]): [description]

    Returns:
        [type]: [description]
    """

    print('Parsing clstr file')
    ### Parse cluster annotation
    proteinclstr = dict()
    bin_proteomes = defaultdict(set)
    motherbin_to_bins = defaultdict(set)

    i = 0
    clstr = None
    cdhit_file = os.path.join(args.v,'eggnog','proteins.nr.faa.clstr')
    with open(cdhit_file,'r') as infile:
        for line in infile:
            if line[0] == '>':
                i += 1 
                clstr = str(i)
            else:
                line = line.strip().split('>')[1]
                proteinid = line.split('...')[0]
                binid = '_'.join( proteinid.split('_')[0:2] )
                motherbin = proteinid.split('_')[1]

                motherbin_to_bins[motherbin].add(binid)
                bin_proteomes[binid].add(proteinid)

                proteinclstr[proteinid] = clstr
    
    return proteinclstr, bin_proteomes, motherbin_to_bins


def load_GO_mapping():
    """[Load GO-slim information for each GO. Plus, get the GO-slim annotation from interpro2go flat-file]

    Returns:
        [type]: [description]
    """

    goslim = {}
    goterm, name, namespace = None, None, None
    dbdir = '/home/projects/cpr_10006/projects/phamb/databases/interpro'
    goslim_files = ['goslim_metagenomics.obo','goslim_agr.obo']
    for gofile in goslim_files:
        with open(os.path.join(dbdir,gofile),'r') as infile:
            for line in infile:
                if line[:3] == 'id:':
                    goterm = line.strip().split(' ')[1]
                    goslim[goterm] = {'name':None,'namespace':None}
                elif line[:5] == 'name:':
                    name = line.strip().split(':')[1]
                    name = name.strip()
                    goslim[goterm]['name'] = name
                elif line[:9] == 'namespace':
                    namespace = line.strip().split(':')[1]
                    namespace = namespace.strip()
                    goslim[goterm]['namespace'] = namespace
                if line[:6] == 'alt_id':
                    altgoterm = line.strip().split(' ')[1]
                    goslim[altgoterm] = {'name':None,'namespace':None}
                    goslim[altgoterm]['name'] = name
                    goslim[altgoterm]['namespace'] = namespace
    goslim = {k:v for k,v in goslim.items() if 'GO' in k}
    
    ### Parse interpro2go annotation
    interpro2go = {}
    with open(os.path.join(dbdir,'interpro2go')) as infile:
        for line in infile:
            if line[0] == '!':
                continue       
            line = line.strip()
            iprid= line.split('>')[0].split(' ')[0]
            iprid = iprid.replace('InterPro:','')
            dircount = line.count('>')
            if dircount == 1:
                high_level_function, GOid = line.split('>')[1].split(';')
                high_level_function = high_level_function.replace('GO:','').strip()
                interpro2go[iprid] = {'GO':GOid, 'GO_low':high_level_function}
            else:
                print('>>>>>ENZYME>>>>>',line)

    return goslim , interpro2go  

def parse_interproscan(args,proteinclstr):
    """[summary]

    Args:
        args ([type]): [description]

    Returns:
        [type]: [description]
    """
    
    interprofile = os.path.join(args.v,'eggnog','interproscan.tsv')

    print('Loading GOslim mapping')
    goslim,interpro2go = load_GO_mapping()

    print('Parsing interprofile')
    interpro_protein = {}
    if not os.path.exists(interprofile):
            print('You need to generate this file first:', interprofile)
            sys.exit(1)
    else:
        parsed_proteins = set()
        with open(interprofile,'r') as infile:
            for rawline in infile:
                line = rawline.strip().split('\t')
                protein = line[0]
                clstr = proteinclstr[protein]
                
                if protein in parsed_proteins:
                    continue
                parsed_proteins.add(protein)

                iprid, iprdesc = None , None 
                if len(line) > 11:
                    iprid = line[11]
                    iprdesc = line[12]
                
                database = line[4]
                databaseid = line[5]
                
                GO_low = None
                if iprid in interpro2go:
                    GO_low = interpro2go[iprid]['GO_low']             
                GO = None
                GO_medium = None
                GO_high = None
                if 'GO:' in rawline:
                    GO = line[13].split('|')                        

                    ### Get higher level annotation for protein
                    for GOid in GO:
                        if GOid in goslim:
                            GO_medium , GO_high = goslim[GOid]['name'], goslim[GOid]['namespace']
                
                interpro_protein[clstr] = {'database':database,'databaseid':databaseid,'GO':GO,'IPR':(iprid,iprdesc),'GO_low':GO_low,'GO_medium':GO_medium,'GO_high':GO_high}
    
    ### Set NA for the remaining protein clstrs 
    for clstr in set(proteinclstr.values()):
        if not clstr in interpro_protein:
            interpro_protein[clstr] = {'database':'None','databaseid':'None','GO':'GO','IPR':('None','None'),'GO_low':'None','GO_medium':'None','GO_high':'None'}
    
    return interpro_protein
    


def parse_DGR_search(args,proteinclstr):
    """[summary]

    Args:
        args ([type]): [description]
        proteinclstr ([type]): [description]

    Returns:
        [type]: [description]
    """

    DGR_file = os.path.join(args.v,'eggnog','DGRsearch.txt')
    DGR_clstr = dict()

    DGR_type = {'Clean_DGR_clade_1':'viral',
    'Clean_DGR_clade_2':'cellular',
    'Clean_DGR_clade_3':'viral/celluar',
    'Clean_DGR_clade_4':'viral',
    'Clean_DGR_clade_5':'cellular',
    'Clean_DGR_clade_6':'viral',
    'RVT_1':'nonDGR',
    'RVT_2':'nonDGR',
    'RVT_3':'nonDGR',
    'RVT_N':'nonDGR'}

    with open(DGR_file,'r') as infile:
        parsed_proteins = set()
        for line in infile:
            if line[0] != '#':
                line = line.strip().split()
                protein, DGR, hmmscore  = line[0], line[2], float(line[5])
                dgrtype =  DGR_type[DGR]

                if hmmscore >= 50 and protein not in parsed_proteins:
                    
                    parsed_proteins.add(protein)
                    clstr = proteinclstr[protein] 
                    
                    if not clstr in DGR_clstr:
                        DGR_clstr[clstr] = (protein, DGR, hmmscore,dgrtype)
                    else:
                        if hmmscore > DGR_clstr[clstr][2]:
                            DGR_clstr[clstr] = (protein, DGR, hmmscore,dgrtype)
    
    return DGR_clstr
                    




def write_bin_protein_annotation(args,binids, annotation_files,fileout):
    """[summary]

    Args:
        args ([type]): [description]
        binids ([type]): [description]
        annotation_files ([type]): [description]
        fileout ([type]): [description]
    """

    interpro_protein = annotation_files['interpro_protein']
    proteinclstr = annotation_files['proteinclstr']
    bin_proteomes = annotation_files['bin_proteomes']
    DGR_clstr = annotation_files['DGR']

    for binid in binids:
        motherbin = binid.split('_')[1]

        for protein in bin_proteomes[binid]:
            nproteins = len(bin_proteomes[binid])
            clstr = proteinclstr[protein]

            DGR = 'NA'
            if clstr in DGR_clstr:
                p, DGR, hmmscore,dgrtype = DGR_clstr[clstr]
                DGR = ':'.join([DGR,dgrtype])

            if clstr in interpro_protein:

            
                interpro_entry = interpro_protein[clstr]
                db = interpro_entry['database']
                dbid = interpro_entry['databaseid']
                iprdesc = interpro_entry['IPR'][1]
                GO = ','.join([ str(interpro_entry['GO_low']), str(interpro_entry['GO_medium']), str(interpro_entry['GO_high']) ]) 

                if db != 'None':
                    lineout = '\t'.join( [protein,str(clstr),binid, motherbin, str(nproteins), db, dbid, str(iprdesc),GO, DGR]  )
                    outfile.write(lineout+'\n')


def write_protein_clusters(args,annotation_files):
    """[summary]

    Args:
        args ([type]): [description]
        annotation_files ([type]): [description]
    """

    protein_clstr_out =  os.path.join(args.o,'proteinclstr.interpros.txt')
    proteinclstr = annotation_files['proteinclstr']
    DGR_clstr = annotation_files['DGR']

    with open(protein_clstr_out,'w') as outfile:
        for clstr in interpro_protein:

            DGR = 'NA'
            if clstr in DGR_clstr:
                p, DGR, hmmscore,dgrtype = DGR_clstr[clstr]
                DGR = ':'.join([DGR,dgrtype])

            interpro_entry = interpro_protein[clstr]
            db = interpro_entry['database']
            dbid = interpro_entry['databaseid']
            iprid = interpro_entry['IPR'][0]
            iprdesc = interpro_entry['IPR'][1]
            GO = ','.join([ str(interpro_entry['GO_low']), str(interpro_entry['GO_medium']), str(interpro_entry['GO_high']) ])
            lineout = '\t'.join( [str(clstr), db, dbid, str(iprid), str(iprdesc),GO, DGR] )
            outfile.write(lineout+'\n')


def get_NC_viral_bins(args):
    """[summary]

    Args:
        args ([type]): [description]
    """
    NC_bins = set()
    MQ_LQ_ND_bins = set()
    with open(os.path.join(args.v,'quality_summary.tsv'),'r') as infile:
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

            quality = line[quality_index]
            genome_copies = float(line[genome_copies_index])
            contaminaton = float(line[contamination_index])
            completeness_method = line[method_index]

            if genome_copies > 1.25 or contaminaton > 10:
                continue

            if quality in ['High-quality','Complete'] and completeness_method =='AAI-based':
                NC_bins.add(binid)
            else:
                MQ_LQ_ND_bins.add(binid)
    
    return NC_bins, MQ_LQ_ND_bins



if __name__ == "__main__":
    args = parser.parse_args()

    proteinclstr, bin_proteomes,motherbin_to_bins  = parse_protein_clstrs(args)
    interpro_protein = parse_interproscan(args, proteinclstr )

    ### Get DGR annotation for each protein-clstr
    DGR_clstr = parse_DGR_search(args,proteinclstr) 

    NC_bins, MQ_LQ_ND_bins = get_NC_viral_bins(args)


    annotation_files = {'interpro_protein':interpro_protein, 
                        'proteinclstr':proteinclstr,
                        'bin_proteomes':bin_proteomes,
                        'DGR':DGR_clstr}

    write_protein_clusters(args,annotation_files)

    if not os.path.exists(args.o):
        os.system('mkdir -p {}'.format(args.o))
    fileout = os.path.join(args.o,'nc.bins.interpros.txt')
    with open(fileout,'w') as outfile:
        outfile.write('protein\tproteinclstr\tbinid\tmotherbin\tproteome\tdb\tdbid\tIPRdesc\tGO\tDGR\n')
        write_bin_protein_annotation(args, NC_bins, annotation_files,outfile)

    fileout = os.path.join(args.o,'MQLQND.bins.interpros.txt')
    with open(fileout,'w') as outfile:
        outfile.write('protein\tproteinclstr\tbinid\tmotherbin\tproteome\tdb\tdbid\tIPRdesc\tGO\tDGR\n')
        write_bin_protein_annotation(args, MQ_LQ_ND_bins, annotation_files,outfile)



