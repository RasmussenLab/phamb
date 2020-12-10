#!/bin/python
import argparse
import os
import sys
import gzip
import taxopy
import json

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-i',help='Directory with spacer-blast files')
parser.add_argument('-d',help='Database of GB accessionids to Taxid ')

def lineage_to_name(LCA_number,taxdb):
        dlineage = {'superkingdom':'NA','kingdom':'NA','phylum':'NA','class':'NA','order':'NA','family':'NA','genus':'NA','species':'NA','subspecies':'NA'}
        current_level = str(LCA_number)
        lineage = taxdb.taxid2name[current_level]
        rank = taxdb.taxid2rank[current_level]
        dlineage[rank] = lineage
        while current_level != '1':
            current_level = taxdb.taxid2parent[current_level]
            lineage  = taxdb.taxid2name[current_level] 
            rank= taxdb.taxid2rank[current_level] 
            dlineage[rank] = lineage

        return dlineage

def read_spacer_dict(args):
    """[Reads a .fna file with spacers and maps the spacer-sequence to the entry naame ]

    Args:
        args ([type]): [description]
    """
    spacerfile = os.path.join(args.i, 'spacers.fna')
    magspacerdict = dict()
    spacerid = None
    with open(spacerfile,'r') as infile:
        for line in infile:
            if line[0] == '>':
                spacerid = line[1:].strip()
            else:
                seq = line.strip()
                magspacerdict[spacerid] = seq
    return magspacerdict




def read_MAG_taxonomy(args):
    """[Read GTDB-TK taxonomy of MAGs]

    Args:
        args ([type]): [description]

    Returns:
        [dict]: [Dictionary with MAG taxonomic lineage]
    """

    #gtdb_file = os.path.join(args.b,'nc_gtdbtk_MQNC/gtdbtk.bac120.summary.tsv')
    gtdb_file = '07_binannotation/bacteria/nc_gtdbtk_MQNC/mag_taxonomy.txt'
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
    return MAG_tax


def get_MAG_overview(args):
    """[Read CheckM file of MQNC MAGs]

    Args:
        args ([type]): [description]
    """
    MAGs = dict()
    checkm_file = '07_binannotation/bacteria/checkm/HMP2.all.MQNC.MAGS'
    with open(checkm_file,'r') as infile:
        header = infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            genomeid = line[0]
            motherbin = genomeid.split('_')[1]
            completeness = float(line[-3])
            qual = None
            if completeness >= 90:
                qual = 'HQ'
            else:
                qual = 'MQ'

            if not motherbin in MAGs:
                MAGs[motherbin] = {}
                MAGs[motherbin][genomeid] = [0,qual]
            else:
                MAGs[motherbin][genomeid] = [0,qual]
    return MAGs 







def load_acc2taxid(args):
    """Load Table with accession to taxid

    Args:
        args ([type]): [description]

    Returns:
        [Dictionary]: [Dict with accession to taxidi]
    """

    acc2taxid = {}
    with gzip.open(args.d,'rt') as infile:
        for line in infile:
            line = line.strip().split()
            acc = line[1]
            taxid = line[2]
            acc2taxid[acc] = taxid
    return acc2taxid


def parse_CRISPR_hits(args):
    """Function for parsing CRISPR-spacer blast-files

    Args:
        args (arguments from commandline): [arguments from commandline]

    """


    ### Parse blastn file of SPACEROME
    spaceromefile = os.path.join(args.i, 'spacerome.viralfragment.m6')
    phage_spacerome = dict()
    with open(spaceromefile,'r') as infile:
        for line in infile:
            line = line.strip().split()

            accession = line[0]
            if '|' in accession:
                accession = accession.split('|')[0]
            
            phagebin = line[1]
            phagemotherbin = phagebin.split('_')[1]
            identity = float(line[2])
            spacer_covered,spacer_len = int(line[3]),int(line[-2])
            bitscore = float(line[-3])
            mismatches = int(line[4])
            if mismatches > 2 or bitscore < 45:
                continue

            if spacer_covered/spacer_len >= 0.6 and identity >= 95:
                if not phagemotherbin in phage_spacerome:
                    phage_spacerome[phagemotherbin] = {}
                    phage_spacerome[phagemotherbin][phagebin] = (accession,identity,bitscore)
                else:
                    if phagebin in phage_spacerome[phagemotherbin]:
                        if bitscore > phage_spacerome[phagemotherbin][phagebin][2]:
                            phage_spacerome[phagemotherbin][phagebin] = (accession,identity,bitscore)
                    else:
                        phage_spacerome[phagemotherbin][phagebin] = (accession,identity,bitscore)
    
    ### Parse MAG spacers
    MAGs = get_MAG_overview(args)
    magspacerdict = read_spacer_dict(args)
    phage_MAG = dict()

    ### Parse blastn file of SPACEROME
    spaceromefile = os.path.join(args.i, 'spacers.viralfragment.m6')
    with open(spaceromefile,'r') as infile:
        for line in infile:
            line = line.strip().split()

            spacername = line[0]
            MAG = '_'.join(spacername.split('_')[:2])
            MAGmotherbin = spacername.split('_')[1]
            phagebin = line[1]
            phagemotherbin = phagebin.split('_')[1]
            identity = float(line[2])
            spacer_covered,spacer_len = int(line[3]),int(line[-2])
            bitscore = float(line[-3])
            mismatches = int(line[4])

            ### Stringent filtering of hits
            if mismatches > 2 or bitscore < 45:
                continue

            if spacer_covered/spacer_len >= 0.6 and identity >= 95:

                ### Update hits dictionary for MAG
                MAGs[MAGmotherbin][MAG][0] = 1

                if not phagemotherbin in phage_MAG:
                    phage_MAG[phagemotherbin] = {}
                    phage_MAG[phagemotherbin][phagebin] = {}
                    phage_MAG[phagemotherbin][phagebin][MAG] = (MAG, identity,bitscore,spacername)
                else:
                    if phagebin in phage_MAG[phagemotherbin]:
                        if MAG in phage_MAG[phagemotherbin][phagebin]:
                            if bitscore > phage_MAG[phagemotherbin][phagebin][MAG][2]:
                                phage_MAG[phagemotherbin][phagebin][MAG] = (MAG, identity,bitscore,spacername)
                        else:
                            phage_MAG[phagemotherbin][phagebin][MAG] = (MAG, identity,bitscore,spacername)
                    else:
                        phage_MAG[phagemotherbin][phagebin] = {}
                        phage_MAG[phagemotherbin][phagebin][MAG] = (MAG, identity,bitscore,spacername)

    MAGtax = read_MAG_taxonomy(args)
    lineage = ['superkingdom','phylum','class','order','family','genus','species']

    ### Write Host taxonomy of each Phage bin 
    tableout = os.path.join(args.i,'phage_host_taxonomy.txt')

    with open(tableout,'w') as outfile:
        for phagemotherbin in phage_MAG:
            for phagebin in phage_MAG[phagemotherbin]:
                for MAG in phage_MAG[phagemotherbin][phagebin]:
                    MAG,identity,bitscore,spacername = phage_MAG[phagemotherbin][phagebin][MAG]
                    seq = magspacerdict[spacername]
                    motherbin = MAG.split('_')[1]
                    taxstring = 'NA'
                    if motherbin in MAGtax:
                        taxstring = ';'.join([ MAGtax[motherbin][lin] for lin in lineage ])

                    outline = [phagemotherbin,phagebin,MAG,str(identity),taxstring,seq,spacername]
                    outfile.write('\t'.join(outline)+'\n')



    ### Write overview with Number of CRISPR annotations to a Phage by Motherbin
    tableout =  os.path.join(args.i,'MAG_overview_taxonomy.txt')
    
    with open(tableout,'w') as outfile:
        for MAG in MAGs:
            nmags = len(MAGs[MAG])
            genomes = list(MAGs[MAG].keys())
            #nspacerhits = 0 
            #for binid in genomes:
            #    if MAGs[MAG][binid] == 1:
            #        nspacerhits += 1
            #perc = round((nspacerhits/nmags)*100,2)
            taxstring = 'NA'
            if MAG in MAGtax:
                    taxstring = ';'.join([ MAGtax[MAG][lin] for lin in lineage ])
            
            for binid in genomes:
                qual = MAGs[MAG][binid][1]
                hit = MAGs[MAG][binid][0]
                outline = [MAG,binid,qual,str(hit),str(nmags),taxstring]
                outfile.write('\t'.join(outline)+'\n')

    ### Write Host taxonomy of each Phage bin 
    print('Loading Accession2taxid')
    acc2taxid = load_acc2taxid(args)
    print('DONE Loading Accession2taxid')

    ### Taxonomy 
    nodes = '/home/projects/cpr_10006/projects/phamb/databases/TaxonomyDB/NCBI_Taxonomy_new_db_2019-09-20/nodes.dmp'
    names ='/home/projects/cpr_10006/projects/phamb/databases/TaxonomyDB/NCBI_Taxonomy_new_db_2019-09-20/names.dmp'
    taxdb = taxopy.TaxDb(nodes_dmp=nodes, names_dmp=names, keep_files=True)
    

    tableout = os.path.join(args.i,'phage_host_spacerome_taxonomy.txt')
    with open(tableout,'w') as outfile:
        for phagemotherbin in phage_spacerome:
            for phagebin in phage_spacerome[phagemotherbin]:
                ACC,identity,bitscore = phage_spacerome[phagemotherbin][phagebin]
                if ACC in acc2taxid:
                    taxid = acc2taxid[ACC]
                    taxonomy = lineage_to_name(taxid,taxdb) 
                    taxstring = ';'.join([ taxonomy[lin] for lin in lineage ])
                else:
                    taxstring = 'NA'
                outline = [phagemotherbin,phagebin,ACC,str(identity),taxstring]
                outfile.write('\t'.join(outline)+'\n')
    



if __name__ == "__main__":
    args = parser.parse_args()

    #accdict = parse_gb_acc2taxid(args)
    parse_CRISPR_hits(args)




