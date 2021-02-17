#!/bin/python
import argparse
import os
import sys
import pathlib
from Bio import SeqIO
from collections import defaultdict
import collections
import re
import numpy as np
import csv


parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-v', help='checkv directory')
parser.add_argument('-i', help='interproscan directory')


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
    cdhit_file = os.path.join(args.v,'virus_proteins.nr.faa.clstr')
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
    
    interprofile = os.path.join(args.i,'virus_proteins.nr.clean.faa.tsv')

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

                if not clstr in interpro_protein:
                    interpro_protein[clstr] = dict()
                
                iprid, iprdesc = None , None 
                if len(line) > 11:
                    iprid = line[11]
                    iprdesc = line[12]
                
                db = line[3]
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
                
                interpro_protein[clstr][db] = {'database':database,'databaseid':databaseid,'GO':GO,'IPR':(iprid,iprdesc),'GO_low':GO_low,'GO_medium':GO_medium,'GO_high':GO_high}
    

    return interpro_protein
    



def parse_DGR_search(args,proteinclstr):
    """[summary]

    Args:
        args ([type]): [description]
        proteinclstr ([type]): [description]

    Returns:
        [type]: [description]
    """

    DGR_file = os.path.join(args.i,'HMP2_DGRsearch.txt')
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

    DGR_clstr = dict()
    DGR_viral_cluster = dict()
    with open(DGR_file,'r') as infile:
        parsed_proteins = set()
        for line in infile:
            if line[0] != '#':
                line = line.strip().split()
                protein, DGR, hmmscore  = line[0], line[2], float(line[5])
                viral_cluster = protein.split('_')[1]
                dgrtype =  DGR_type[DGR]

                ### Cutoff suggested in S. Roux article 2020
                if hmmscore >= 30:                    
                    parsed_proteins.add(protein)
                    
                    if not protein in DGR_clstr:
                        DGR_clstr[protein] = (protein, DGR, hmmscore,dgrtype)
                    else:
                        if hmmscore > DGR_clstr[protein][2]:
                            DGR_clstr[protein] = (protein, DGR, hmmscore,dgrtype)
                    
                    ### 
                    if not viral_cluster in DGR_viral_cluster:
                        DGR_viral_cluster[viral_cluster] = dict()
                        DGR_viral_cluster[viral_cluster][protein] = (DGR,hmmscore)
                    else:

                        if protein in DGR_viral_cluster[viral_cluster]:
                            if hmmscore > DGR_viral_cluster[viral_cluster][protein][1]:
                                DGR_viral_cluster[viral_cluster][protein] = (DGR,hmmscore)
                        else:
                            DGR_viral_cluster[viral_cluster][protein] = (DGR,hmmscore)

    DGR_fileout = os.path.join(args.i, 'DGR.RT.scan.tsv')
    with open(DGR_fileout,'w') as out:
        for viral_cluster in DGR_viral_cluster:
            for protein in DGR_viral_cluster[viral_cluster]:
                DGRtype, hmmscore = DGR_viral_cluster[viral_cluster][protein]
                lineout = [protein, viral_cluster, DGRtype, DGR_type[DGRtype], hmmscore]
                lineout = [str(x) for x in lineout]
                out.write('\t'.join(lineout)+'\n')


    
    DGR_search_directory = '/home/projects/cpr_10006/projects/phamb/databases/dgr_scripts/DGR_identification/HMP2/virus'
    DGR_target_fileout = os.path.join(args.i, 'DGR.Target.scan.tsv')
    
    DGR_target_scan = os.path.join(DGR_search_directory,'interproscan','target_sequences_clean.faa.tsv')
    DGR_target_PC = os.path.join(DGR_search_directory,'target_clustering','All_RT_1k_DGR_target_detection_filtered_targets_prot_to_PC.tsv')

    class PC_class:
        def __init__(self):
            self.proteins = []
            self.annotations = {'Pfam':[],'Gene3D':[],'TIGRFAM':[],'SUPERFAMILY':[]}
            self.majority_annotation = dict()
            self.nproteins = 0
            self.DGRtype = {k:0 for k in list(DGR_type)}

    protein_to_PC = dict()
    PCs = dict()
    with open(DGR_target_PC,'r') as infile:
        for line in infile:
            PCid, protein = line.strip().split('\t')
            protein_to_PC[protein] = PCid
            if not PCid in PCs:
                PC = PC_class()
                PC.proteins += [protein]
                PC.nproteins += 1
                PCs[PCid] = PC
            else:
                PCs[PCid].proteins += [protein]
                PC.nproteins += 1
    
    ### Parse Interpro     
    interpro_protein = dict()
    with open(DGR_target_scan,'r') as infile:
        for rawline in infile:
            line = rawline.strip().split('\t')
            protein = line[0]
            viral_cluster = protein.split('_')[1] 
            if not protein in protein_to_PC:
                continue
            PCid = protein_to_PC[protein]
            if not PCid in interpro_protein:
                interpro_protein[PCid] = dict()
            iprid, iprdesc = None , None 
            if len(line) > 11:
                iprid = line[11]
                iprdesc = line[12]
            db = line[3]
            databaseid = line[4]
            databasedesc = line[5]
            PCs[PCid].annotations[db] += [(protein,databaseid,databasedesc,iprid,iprdesc)]

    ### Calculate majority vote of PC interproscan annotations
    for PCid in PCs:
        for db in PCs[PCid].annotations:
            entries = PCs[PCid].annotations[db]
            nentries = len(entries)
            dbids = [ x[1] for x in entries]
            counts = collections.Counter(dbids)
            counts_frequencies = {k: round(v/nentries,2) for k,v in counts.items()}
            if len(counts_frequencies) != 0:
                dominant_annotation = [k for k, v in sorted(counts_frequencies.items(), key=lambda item: item[1],reverse=True)][0]
                PCs[PCid].majority_annotation[db] = (dominant_annotation, counts_frequencies[dominant_annotation],nentries )

            

    with open(DGR_target_fileout,'w') as outfile:
        outfile.write('protein\tproteinclstr\tbinid\tmotherbin\tinterprodb\tdb\tdbid\tIPRdesc\n')
        for protein in dgr_clstrs:
            motherbin = protein.split('_')[1]
            binid = '_'.join(protein.split('_')[:2])
            clstr = dgr_clstrs[protein]
            if clstr in interpro_protein:
                for interprodb in interpro_protein[clstr]: 
                    interpro_entry = interpro_protein[clstr][interprodb]
                    db = interpro_entry['database']
                    dbid = interpro_entry['databaseid']
                    iprdesc = interpro_entry['IPR'][1]

                    if db != 'None':
                        lineout = '\t'.join( [protein,str(clstr),binid, motherbin, interprodb,db, dbid, str(iprdesc)]  )
                        outfile.write(lineout+'\n')
    
    ### Depedencies for Protein clustering
    # module load hhsuite/3.3.0  (module load anaconda3/4.4.0 perl/5.24.0 openmpi/gcc/64/1.10.2 )
    # Hdd to change --cpu 1 to --cpu2 in hhsuitedb.py in S. Rouxs script



def write_protein_clusters(args,annotation_files):
    """[summary]

    Args:
        args ([type]): [description]
        annotation_files ([type]): [description]
    """

    protein_clstr_out =  os.path.join(args.v,'proteinclstr.interpros.txt')
    proteinclstr = annotation_files['proteinclstr']
    DGR_clstr = annotation_files['DGR']

    with open(protein_clstr_out,'w') as outfile:
        for clstr in interpro_protein:

            if clstr in DGR_clstr:
                p, DGR, hmmscore,dgrtype = DGR_clstr[clstr]
                DGR = ':'.join([DGR,dgrtype])

            interpro_entry = interpro_protein[clstr]
            db = interpro_entry['database']
            dbid = interpro_entry['databaseid']
            iprid = interpro_entry['IPR'][0]
            iprdesc = interpro_entry['IPR'][1]
            GO = ','.join([ str(interpro_entry['GO_low']), str(interpro_entry['GO_medium']), str(interpro_entry['GO_high']) ])
            lineout = '\t'.join( [str(clstr), db, dbid, str(iprid), str(iprdesc),GO] )
            outfile.write(lineout+'\n')


def get_viral_bins(args):
    """[summary]

    Args:
        args ([type]): [description]
    """
    #C_bins = set()
    #MQ_LQ_ND_bins = set()

    quality_file = os.path.join(args.v,'quality_summary.tsv')
    viruses = dict()
    with open(quality_file,'r') as infile:
        reader = csv.DictReader(infile,delimiter='\t') 
        for row in reader:
            contigid  = row['contig_id']
            if not row['checkv_quality'] in ['Medium-quality','High-quality','Complete']:
                continue
            viruses[contigid] = row
    
    return viruses

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

    for binid in binids:
        motherbin = binid.split('_')[1]
        for protein in bin_proteomes[binid]:
            nproteins = len(bin_proteomes[binid])
            clstr = proteinclstr[protein]
            if clstr in interpro_protein:                    
                for interprodb in interpro_protein[clstr]: 
                    interpro_entry = interpro_protein[clstr][interprodb]
                    db = interpro_entry['database']
                    dbid = interpro_entry['databaseid']
                    iprdesc = interpro_entry['IPR'][1]
                    GO = ','.join([ str(interpro_entry['GO_low']), str(interpro_entry['GO_medium']), str(interpro_entry['GO_high']) ]) 

                    if db != 'None':
                        lineout = '\t'.join( [protein,str(clstr),binid, motherbin, str(nproteins), interprodb,db, dbid, str(iprdesc),GO]  )
                        outfile.write(lineout+'\n')
            else:
                lineout = '\t'.join( [protein,str(clstr),binid, motherbin, str(nproteins), 'NA','NA', 'NA', 'NA','NA']  )
                outfile.write(lineout+'\n')



class meh:
    def __init__(self,v,i):
        self.v = v
        self.i = i 

args = meh('07_binannotation/checkv/VAMB_bins','07_binannotation/checkv/interproscan')


if __name__ == "__main__":
    args = parser.parse_args()

    proteinclstr, bin_proteomes,motherbin_to_bins  = parse_protein_clstrs(args)
    interpro_protein = parse_interproscan(args, proteinclstr )

    ### Get DGR annotation for each protein-clstr
    #DGR_clstr = parse_DGR_search(args,proteinclstr) 

    viruses = get_viral_bins(args)


    annotation_files = {'interpro_protein':interpro_protein, 
                        'proteinclstr':proteinclstr,
                        'bin_proteomes':bin_proteomes}

    #write_protein_clusters(args,annotation_files)

    if not os.path.exists(args.i):
        os.system('mkdir -p {}'.format(args.i))
    fileout = os.path.join(args.i,'virus.interpros.txt')
    with open(fileout,'w') as outfile:
        outfile.write('protein\tproteinclstr\tbinid\tmotherbin\tproteome\tinterprodb\tdb\tdbid\tIPRdesc\tGO\n')
        write_bin_protein_annotation(args, viruses, annotation_files,outfile)

    #fileout = os.path.join(args.o,'MQLQND.bins.interpros.txt')
    #with open(fileout,'w') as outfile:
     #   outfile.write('protein\tproteinclstr\tbinid\tmotherbin\tproteome\tdb\tdbid\tIPRdesc\tGO\tDGR\n')
      #  write_bin_protein_annotation(args, MQ_LQ_ND_bins, annotation_files,outfile)



