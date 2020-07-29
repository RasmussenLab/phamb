#!/usr/bin/env python3

import argparse
import os
import sys
import pathlib
import taxopy
from Bio import SeqIO
import subprocess
import copy
import gzip

### custom modulees
import tax_annotations
import checkv_parsers

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-c',help='VAMB clusterfile/clusters.tsv')
parser.add_argument('-a',help='Annotation summary files directory')
parser.add_argument('-v', help='checkv directory')


class motherbintax:
    def __init__(self,motherbinname):
        self.motherbinname = motherbinname
        self.bins = set()
        self.contigs = set()
        self.taxlevels = {0:[],1:[],2:[],3:[],4:[],5:[],6:[]}
        self.taxlevels_freq = {0:None,1:None,2:None,3:None,4:None,5:None,6:None}
        self.taxlevels_LCA = {0:None,1:None,2:None,3:None,4:None,5:None,6:None}
        self.motherbintax = None
        self.motherbintax_GCA = None
        self.motherbintax_VOGclade = None

    def summarise_tax_frequencies(self):
        for i in self.taxlevels.keys():
            
            item_freq = []
            items = self.taxlevels[i]
            non_NA_items = [x for x in items if x != 'NA' ]
            n_NA_items = items.count('NA')

            if len(items) == 0:
                self.taxlevels_LCA[i] =  ('NA','NA')
                continue

            ### If 90% of all hits are NA - let it be NA .
            if n_NA_items/len(items) >= 0.9:
                sorted_item_freq = [ ('NA',float(1)) ] 
                self.taxlevels_LCA[i] = ('NA','NA')
            else:
                unique_elements = set(non_NA_items)
                n_non_NA_items = len(non_NA_items)
                for y in unique_elements:
                    freq = items.count(y)/n_non_NA_items
                    item_freq.append((y,freq))
            
                sorted_item_freq = sorted(item_freq, key = lambda x: x[1], reverse=True)
                self.taxlevels_freq[i] = sorted_item_freq

            ### Determine the LCA - what the highest level of taxonomic consensus
            best_tax = None
 
            if len(sorted_item_freq) == 1:
                best_tax = sorted_item_freq[0]
            else:
                tup = sorted_item_freq
                
                ### Either choose the single taxonomy in the level 
                ### Or choose the taxonomy annotated X% of the contigs based on the non-NA items.
                if len(tup) == 1:
                    best_tax = tup[0]
                elif tup[0][1] >= 0.50:
                    best_tax = tup[0]
                else:
                    best_tax = ('NA','NA')
            self.taxlevels_LCA[i] = best_tax
        
        ### Make LCA consensus string
        taxs = []
        for i in self.taxlevels_LCA.keys():
            leveltax = self.taxlevels_LCA[i][0]
            taxs.append(leveltax)
            self.motherbintax = ';'.join(taxs)
  

class bintax:
    def __init__(self,bin_name):
        self.bin_name = bin_name
        self.contigs = set()
        self.contig_lengths = []
        self.taxlevels = {0:[],1:[],2:[],3:[],4:[],5:[],6:[]}
        self.taxlevels_freq = {0:None,1:None,2:None,3:None,4:None,5:None,6:None}
        self.taxlevels_LCA = {0:None,1:None,2:None,3:None,4:None,5:None,6:None}
        self.bintax = None
        self.GCAtax = None
        self.VOGclade = None
        self.checkv_results = {'completeness':None,'prophage':None,'checkv_quality':None,'miuvig_quality':None,'completeness_method':None,'HMMname':None}
        self.checkv_reference = None
        self.checkv_AAI_confidence = None
     
    def summarise_tax_frequencies(self):
        for i in self.taxlevels.keys():
            
            item_freq = []
            items = self.taxlevels[i]
            non_NA_items = [x for x in items if x != 'NA' ]
            n_NA_items = items.count('NA')

            if len(items) == 0:
                self.taxlevels_LCA[i] =  ('NA','NA')
                continue

            ### If 90% of all hits are NA - let it be NA .
            if n_NA_items/len(items) >= 0.9:
                sorted_item_freq = [ ('NA',float(1)) ] 
                self.taxlevels_LCA[i] = ('NA','NA')
            else:
                unique_elements = set(non_NA_items)
                n_non_NA_items = len(non_NA_items)
                for y in unique_elements:
                    freq = items.count(y)/n_non_NA_items
                    item_freq.append((y,freq))
            
                sorted_item_freq = sorted(item_freq, key = lambda x: x[1], reverse=True)
                self.taxlevels_freq[i] = sorted_item_freq

            ### Determine the LCA - what the highest level of taxonomic consensus
            best_tax = None

            
            if len(sorted_item_freq) == 1:
                best_tax = sorted_item_freq[0]
            else:
                tup = sorted_item_freq
                
                ### Either choose the single taxonomy in the level 
                ### Or choose the taxonomy annotated X% of the contigs based on the non-NA items.
                if len(tup) == 1:
                    best_tax = tup[0]
                elif tup[0][1] >= 0.50:
                    best_tax = tup[0]
                else:
                    best_tax = ('NA','NA')
            self.taxlevels_LCA[i] = best_tax
        
        ### Make LCA consensus string
        taxs = []
        for i in self.taxlevels_LCA.keys():
            leveltax = self.taxlevels_LCA[i][0]
            taxs.append(leveltax)
            self.bintax = ';'.join(taxs)



### Taxonomy 
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

nodes = '/home/projects/cpr_10006/projects/phamb/databases/TaxonomyDB/NCBI_Taxonomy_new_db_2019-09-20/nodes.dmp'
names ='/home/projects/cpr_10006/projects/phamb/databases/TaxonomyDB/NCBI_Taxonomy_new_db_2019-09-20/names.dmp'
taxdb = taxopy.TaxDb(nodes_dmp=nodes, names_dmp=names, keep_files=True)
name2taxid = {v:k for k,v in taxdb.taxid2name.items()}

### Related to checkV files

def parse_checkv_write_overview(args):

      
    def parse_clusterfile(args):
        '''Return list of contigids, size and affiliation'''

        contigs = {}
        clusterfile = args.c
        motherbins = {}
        bins = {}
        with open(clusterfile,'r') as infile:
            for line in infile:
                cluster, contig = line.strip().split()
                contigs[contig] = {'VOGenrichment':None,'contiglength':None}
                sample = contig.split('_')[0]
                binid = sample + '_' + cluster
                
                ### make object for Bin 
                if binid not in bins:
                    binobject = bintax(binid)
                else:
                    binobject = bins[binid]

                ### make object for cluster/motherbin
                if cluster not in motherbins:
                    motherbin = motherbintax(cluster)
                else:
                    motherbin = motherbins[cluster]

                motherbin.contigs.add(contig)
                motherbin.bins.add(binid)

                binobject.contigs.add(contig)
                length = contig.split('length')[1].split('_')[1]
                contigs[contig]['contiglength'] = float(length)
                binobject.contig_lengths += [float(length)]
                bins[binid] = binobject
                motherbins[cluster] = motherbin

        return bins, contigs, motherbins


    def parse_VOGs(contigs):
        '''
        Calculate Tax consensus for each contig if it have more than 2 different VOG hits
        '''

        def load_VOG_table():
            taxlevels = ['superkingdom','phylum','class','order','family','genus','species']
            vogtax = '/home/projects/cpr_10006/projects/phamb/databases/VOG/vog.lca.tax.tsv' 
            vog_tax = dict()

            if not os.path.exists(vogtax):
                print('Missing VOG database file !')
                sys.exit(1)

            with open(vogtax, 'r') as infile:
                for line in infile:
                    line = line.strip().split('\t')
                    vog = line[0]
                    vog_tax[vog] = {i:None for i in range(len(taxlevels))}
                    for i in range(len(line[1:])):
                        vog_tax[vog][i] = line[i+1]
            return vog_tax


        def vog_consensus(vogs,vog_tax):
            taxlevels = {0:[],1:[],2:[],3:[],4:[],5:[],6:[]}
            for vog in vogs:
                taxinfo = vog_tax[vog]
                for i in range(7):
                    taxlevelentry = taxinfo[i]
                    taxlevels[i] += [taxlevelentry]
            #print(taxlevels)
            taxlevels_freq = {0:None,1:None,2:None,3:None,4:None,5:None,6:None}
            taxlevels_LCA = {0:None,1:None,2:None,3:None,4:None,5:None,6:None}

            for i in taxlevels.keys():
                item_freq = []
                items = taxlevels[i]
                non_NA_items = [x for x in items if x != 'NA' ]
                n_NA_items = items.count('NA')
                
                if len(items) == 0:
                    taxlevels_LCA[i] =  ('NA','NA')
                    continue
                if n_NA_items/len(items) >= 0.9:
                    sorted_item_freq = [ ('NA',float(1)) ] 
                    taxlevels_LCA[i] = ('NA','NA')
                else:
                    unique_elements = set(non_NA_items)
                    n_non_NA_items = len(non_NA_items)
                    for y in unique_elements:
                        freq = items.count(y)/n_non_NA_items
                        item_freq.append((y,freq))
                
                    sorted_item_freq = sorted(item_freq, key = lambda x: x[1], reverse=True)
                    taxlevels_freq[i] = sorted_item_freq

                ### Determine the LCA - what the highest level of taxonomic consensus
                best_tax = None

                if len(sorted_item_freq) == 1:
                    best_tax = sorted_item_freq[0]
                else:
                    tup = sorted_item_freq
                        
                    ### Either choose the single taxonomy in the level 
                    if len(tup) == 1:
                        best_tax = tup[0]
                    elif tup[0][1] >= 0.5:
                        best_tax = tup[0]
                    else:
                        best_tax = ('NA','NA')
                taxlevels_LCA[i] = best_tax

            return taxlevels_LCA

        vog_tax = load_VOG_table()
        vogfile = os.path.join(args.a,'hmmVOG.allsamples.txt')
        if not os.path.exists(vogfile):
            print('Missing VOG file !')
            sys.exit(1)
        contig_vog_tax = {}
        vogcontig = {}
        with open(vogfile,'r') as infile:
            for line in infile:
                contig, protein, vog, score, evalue = line.strip().split('\t')
                if contig not in contigs:
                    continue
                if contig not in vogcontig:
                    vogcontig[contig] = []
                    vogcontig[contig] += [vog]
                else:
                    vogcontig[contig] += [vog]
        
        for contig in vogcontig:
            contiglength = contigs[contig]['contiglength']
            vogs = vogcontig[contig]
            nvogs = len(vogs)
            vogs_per_kb = nvogs/(contiglength/1000)
            if nvogs >= 2:
                contigs[contig]['VOGenrichment'] = round(vogs_per_kb,2)

                ### assess contig taxonomy based on vogs 
                contig_vog_consensus = vog_consensus(vogs,vog_tax)
                contig_vog_tax[contig] = contig_vog_consensus


        return contigs, contig_vog_tax
    
        
    def bin_VOG_tax(bins,motherbins,contig_vog_tax):
        '''
        Based on the contigs in each Bin - try to assign Taxonomy using VOG HMM hits 
        '''
        
        ### SUmmarise tax pr. Bin 
        for binid in bins:
            binobject = bins[binid]
            for contig in binobject.contigs:
                if contig in contig_vog_tax:
                    for i in range(7):
                        taxentry = contig_vog_tax[contig][i][0]
                        binobject.taxlevels[i] += [taxentry]
            binobject.summarise_tax_frequencies()
            bins[binid] = binobject

        ### Do the same for all contigs belonging to the cluster/motherbin
        for cluster in motherbins:
            binobject = motherbins[cluster]
            for contig in binobject.contigs:
                if contig in contig_vog_tax:
                    for i in range(7):
                        taxentry = contig_vog_tax[contig][i][0]
                        binobject.taxlevels[i] += [taxentry]
            binobject.summarise_tax_frequencies()
            motherbins[cluster] = binobject
        
        return bins, motherbins

    

    ### stuff to do
    bins, contigs, motherbins = parse_clusterfile(args)

    ### VOG enrichment analysis
    
    contigs, contig_vog_tax = parse_VOGs(contigs)

    ### Summarise VOG Taxonomy for each Bin
    bins,motherbins = bin_VOG_tax(bins,motherbins,contig_vog_tax)

    del contigs
    del contig_vog_tax 

    ### Load checkV results
    checkv_db_dir = '/home/projects/cpr_10006/projects/phamb/databases/checkv/checkv-db-v0.6/genome_db'
    checkv_genome_gcatax, gcatax, DTRtax  = checkv_parsers.load_checkv_db_files(checkv_db_dir,name2taxid,taxdb)
    bins = checkv_parsers.read_checkv_quality(args,bins)
    bins,motherbins = checkv_parsers.read_checkv_completeness(args,bins,motherbins,checkv_genome_gcatax,gcatax,DTRtax)

    ### Write final overview of Viral Genomes
    checkv_parsers.write_bin_overview(args,bins,motherbins)


### Run annotation programmes on NC Viral Genomes
def annotation_ncgenomes_proteomes(args):
    '''
    Remember to load the following modules 
    module load perl/5.24.0 ncbi-blast/2.8.1+ hmmer/3.2.1 diamond/0.9.29
    '''

    
    checkv_directory = args.v
    prodigalfile_AA = os.path.join(checkv_directory,'nc_genomes_proteins.faa')
    
    
    ### BlastP against Crassphage
    blastp_crass = os.path.join(checkv_directory,'nc_genomes.crassphage.polterm.m6')
    databasedirectory = '/home/projects/cpr_10006/projects/phamb/databases/'
    blastp_executable= '/services/tools/ncbi-blast/2.8.1+/bin/blastp'
    if os.path.exists(blastp_crass):
        print('BlastP file allready created:',blastp_crass,' delete file to recreate it')
    else:
        tax_annotations.crass_blast(databasedirectory, blastp_executable, prodigalfile_AA, blastp_crass)


    ### BlastP against Eukaryotic Viruses 
    databasedirectory = '/home/projects/cpr_10006/projects/phamb/databases/'
    diamond_executable = '/services/tools/diamond/0.9.29/diamond'
    blastp_euk = os.path.join(checkv_directory,'nc_genomes.eukaryotic.viruses.m6')
    genomes = os.path.join(checkv_directory,'nc_genomes.fna')
    if os.path.exists(blastp_euk):
        print('BlastP file allready created:',blastp_euk,' delete file to recreate it')
    else:
        tax_annotations.herpes_blast(databasedirectory,diamond_executable, genomes, blastp_euk)


    ### BlastP ALL contigs against Eukaryotic Viruses -
    blastp_euk_all = os.path.join(checkv_directory,'cleaned_genomes.eukaryotic.viruses.m6')
    genomes = os.path.join(checkv_directory,'cleaned_contigs.fna')
    if os.path.exists(blastp_euk_all):
        print('BlastP file allready created:',blastp_euk_all,' delete file to recreate it')
    else:
        tax_annotations.herpes_blast(databasedirectory,diamond_executable, genomes, blastp_euk_all)

    ### Search against RVDB HMM profiles
    databasedirectory = '/home/projects/cpr_10006/projects/phamb/databases/'
    HMM_executable= '/services/tools/hmmer/3.2.1/bin/hmmsearch'
    rvdb_ref = os.path.join(checkv_directory,'nc_genomes.RVDB.hmmsearch')
    if os.path.exists(rvdb_ref):
        print('Hmmsearch file allready created:',rvdb_ref,' delete file to recreate it')
    else:
        tax_annotations.RVDB_search(databasedirectory,HMM_executable,prodigalfile_AA, rvdb_ref)


    ### Search against VOG HMM profiles
    # Based on this search - it will be much easier to write out Hallmark proteins for Phylogenetic Trees
    VOG_ref = os.path.join(checkv_directory,'nc_genomes.VOG.hmmsearch')
    if os.path.exists(VOG_ref):
        print('Hmmsearch file allready created:',VOG_ref,' delete file to recreate it')
    else:
        tax_annotations.VOG_search(databasedirectory,HMM_executable, prodigalfile_AA, VOG_ref)


    ### Usearch against UNIprot Viral Proteins
    usearch_executable = '/services/tools/usearch/11.0.667/bin/usearch'
    usearch_trembl = os.path.join(checkv_directory,'nc_genomes.trembl_viral.txt')
    if os.path.exists(usearch_trembl):
        print('Usearch file allready created:',usearch_trembl,' delete file to recreate it')
    else:
        tax_annotations.uniprot_search(databasedirectory,usearch_executable, prodigalfile_AA, usearch_trembl)


### Module of Functions to write an additional Taxonomy annotation for NC viral genomes
def write_extended_taxonomy(args):

    class uniprotbintax:
        def __init__(self,bin_name):
            self.bin_name = bin_name
            self.taxlevels = {0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[]}
            self.taxlevels_freq = {0:None,1:None,2:None,3:None,4:None,5:None,6:None,7:None}
            self.taxlevels_LCA = {0:None,1:None,2:None,3:None,4:None,5:None,6:None,7:None}
            self.VOGclade = None
            self.completeness = None
            self.clustertax = None
            self.uniprottax = None
        
        def summarise_tax_frequencies(self):
            for i in self.taxlevels.keys():
                
                item_freq = []
                items = self.taxlevels[i]
                non_NA_items = [x for x in items if x != 'NA' ]
                n_NA_items = items.count('NA')

                if len(items) == 0:
                    self.taxlevels_LCA[i] =  ('NA','NA')
                    continue

                ### If 90% of all hits are NA - let it be NA .
                if n_NA_items/len(items) >= 0.9:
                    sorted_item_freq = [ ('NA',float(1)) ] 
                    self.taxlevels_LCA[i] = ('NA','NA')
                else:
                    unique_elements = set(non_NA_items)
                    n_non_NA_items = len(non_NA_items)
                    for y in unique_elements:
                        freq = items.count(y)/n_non_NA_items
                        item_freq.append((y,freq))
                
                    sorted_item_freq = sorted(item_freq, key = lambda x: x[1], reverse=True)
                    self.taxlevels_freq[i] = sorted_item_freq

                ### Determine the LCA - what the highest level of taxonomic consensus
                best_tax = None

                
                if len(sorted_item_freq) == 1:
                    best_tax = sorted_item_freq[0]
                else:
                    tup = sorted_item_freq
                    
                    ### Either choose the single taxonomy in the level 
                    ### Or choose the taxonomy annotated X% of the contigs based on the non-NA items.
                    if len(tup) == 1:
                        best_tax = tup[0]
                    elif tup[0][1] >= 0.35:
                        best_tax = tup[0]
                    else:
                        best_tax = ('NA','NA')
                self.taxlevels_LCA[i] = best_tax
            
            ### Make LCA consensus string
            taxs = []
            for i in self.taxlevels_LCA.keys():
                leveltax = self.taxlevels_LCA[i][0]
                taxs.append(leveltax)
                self.uniprottax = ';'.join(taxs)


    def read_rvdb_lca(args):
        '''

        
        Familyid -> LCA 
        LCA -> Lineage 
        Lineage -> Taxonomy  {    Use taxdb for this! Great tool <:o)   }

        Read in LCA of all Families of RVDB database
        structure is 
        annotation/6.txt
        FAM000006
        '''

        rvdb_db =  os.path.join('/home/projects/cpr_10006/projects/phamb/databases/','rvdb/rvdb_family_lca.txt')
        
        ### Make FamilyID to LCA dict
        FAMLCA = dict()
        with open(rvdb_db,'r') as infile:
            infile.readline()
            for line in infile:
                line = line.strip().split()
                famid = int(line[1])
                FAMLCA[famid] = line[4]
                            
        return FAMLCA
    
                
    def read_uniprot_db(args,taxdb):
        databasedirectory = '/home/projects/cpr_10006/projects/phamb/databases/'

        tablefile = os.path.join(databasedirectory,'Demovir/UniprotId_to_taxa.txt')
        #tablefile = '/home/projects/cpr_10006/projects/phamb/databases/Demovir/UniprotId_to_taxa.txt'

        name2taxid = {v:k for k,v in taxdb.taxid2name.items()}

        uniprot_taxonomy = {}
        with open(tablefile,'r') as infile:
            infile.readline()
            for line in infile:
                row, uniprotid, order, family, species = line.strip().split('\t')
                if uniprotid not in uniprot_taxonomy:
                    lowesttaxid = None
                    if species in name2taxid:
                        lowesttaxid = ('species',name2taxid[species])
                    elif family in name2taxid:
                        lowesttaxid = ('family',name2taxid[family])
                    elif order in name2taxid:
                        lowesttaxid = ('order',name2taxid[order])
                    taxid = lowesttaxid[1]
                    taxonomy = lineage_to_name(taxid,taxdb)                    
                    uniprot_taxonomy[uniprotid] = taxonomy

        return uniprot_taxonomy

    def read_crassphage_hits(args):
        '''
        Inspired from "Biology and Taxonomy of crAss-like Bacteriophages, the Most Abundant Virus in the Human Gut"
        '''    
        crassblastfile = os.path.join(args.v,'nc_genomes.crassphage.polterm.m6')
        crassbins = {}
        crassproteins = {'YP_009052554.1':'Terminase_large_subunit','YP_009052497.1':'DNA_polB'}
        with open(crassblastfile,'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                proteinid = line[0]
                binid = '_'.join(proteinid.split('_')[:2])
                proteinhit = crassproteins[line[1]]
                seqid = float(line[2])
                aln = int(line[3])
                evalue = float(line[10])
                if aln < 350 or evalue > 1E-05:
                    continue
                if binid not in crassbins:
                    crassbins[binid] = []
                    crassbins[binid] += [proteinhit,seqid,aln]
                else:
                    crassbins[binid] += [proteinhit,seqid,aln]
        return crassbins

    def read_checkv_quality(args):
        
        checkvfile = os.path.join(args.v,'checkv_taxonomy.txt')
        bins = {}
        vogclades = set()
        with open(checkvfile,'r') as infile:
            infile.readline()
            for line in infile:
                line = line.strip().split()
                binid = line[0]
                quality = line[2]

                if quality =='Complete':
                    completeness = '100'
                else:
                    completeness = line[3]
                if quality =='Not-determined':
                    continue           
                clustertaxonomy = line[6]
                VOGclade = line[7]

                if float(completeness) >= 80:
                    binobject = uniprotbintax(binid)
                    binobject.clustertax = clustertaxonomy
                    binobject.completeness = completeness
                    binobject.VOGclade = VOGclade
                    vogclades.add(VOGclade)
                    bins[binid] = binobject
        return bins, vogclades

    def parse_uniprot_hits(args,taxdb,bins):

        ### Load uniprot db
        uniprot_taxonomy = read_uniprot_db(args,taxdb)

        lineage = ['superkingdom','phylum','class','order','family','genus','species','subspecies']
        uniprot_hits = os.path.join(args.v,'nc_genomes.trembl_viral.txt')
        uniprot_bins = copy.deepcopy(bins)
        if os.path.exists(uniprot_hits):
            with open(uniprot_hits,'r') as infile:
                for line in infile:
                    line = line.strip().split('\t')
                    protid = line[0].split(' ')[0]
                    binid = '_'.join(protid.split('_')[:2])
                    if not binid in uniprot_bins:
                        continue
                    binobject = uniprot_bins[binid]
                    uniprotid = line[1]
                    taxonomy = uniprot_taxonomy[uniprotid]

                    for i, lin in enumerate(lineage):
                        binobject.taxlevels[i] += [ taxonomy[lin] ]
                    uniprot_bins[binid] = binobject
        else:
            print(uniprot_hits,' is missing?')

        uniprot_tax = {}
        for binid in uniprot_bins:
            binobject = uniprot_bins[binid]
            binobject.summarise_tax_frequencies()
            tax = binobject.uniprottax
            uniprot_tax[binid] = tax
        return uniprot_tax


    def parse_rvdb_hits(args,taxdb,bins):

        lineage = ['superkingdom','phylum','class','order','family','genus','species','subspecies']

        FAMLCA = read_rvdb_lca(args)
        rvdb_hmmfile = os.path.join(args.v,'nc_genomes.RVDB.hmmsearch')
        rvdb_bins = copy.deepcopy(bins)
        if os.path.exists(rvdb_hmmfile):
            with open(rvdb_hmmfile,'r') as infile:
                for line in infile:
                    if line[0] == '#':
                        continue
                    line = line.strip().split()
                    score = line[5] 
                    if float(score) >= 50:       
                        proteinid = line[0]
                        binid  = '_'.join( proteinid.split('_')[:2])

                        if not binid in rvdb_bins:
                            continue

                        binobject = rvdb_bins[binid]
                        famid = line[2]
                        famid = int(famid.replace('FAM',''))
                        EVAL = line[4]
                        LCA = int(FAMLCA[famid])
                        if LCA == 0:
                            continue

                        if not LCA in taxdb.taxid2name:
                            continue

                        taxonomy = lineage_to_name(LCA,taxdb)

                        for i, lin in enumerate(lineage):
                            binobject.taxlevels[i] += [ taxonomy[lin] ]
                        rvdb_bins[binid] = binobject
        else:
            print(rvdb_hmmfile,' is missing?')

        rvdb_tax = {}
        for binid in rvdb_bins:
            binobject = rvdb_bins[binid]
            binobject.summarise_tax_frequencies()
            tax = binobject.uniprottax
            rvdb_tax[binid] = tax
        return rvdb_tax
                

    def assign_nr_taxonomy(args):
        
        def load_gi2_name_eukvirus(args):
            gi2name = {}
            viral_euk_proteome = '/home/projects/cpr_10006/projects/phamb/databases/herpes_db/proteomes/herpes.faa'
            with open(viral_euk_proteome,'r') as infile:
                for line in infile:
                    if line[0] == '>':
                        line = line[1:].strip().split('[')
                        gi = line[0].split(' ')[0]
                        if len(line) >= 2: 
                            taxname = line[1].replace(']','')
                            gi2name[gi] = taxname
            
            gi2taxid = {k:name2taxid[v] for k,v in gi2name.items() if v in name2taxid}
            return gi2taxid


        gi2taxid = load_gi2_name_eukvirus(args)

            
        blastp_file = '/home/projects/cpr_10006/projects/phamb/HMP2/07_binannotation/checkv/bins_checkv_out/cleaned_genomes.eukaryotic.viruses.m6'
        bintaxonomy = dict()
        skips = 0
        with open(blastp_file,'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                gi = line[1]
                if gi not in gi2taxid:
                    skips += 1
                    continue
                binid = line[0]
                protein_sim = float(line[2])
                evalue = float(line[-2])

                if protein_sim >= 50 and evalue <= 1e-05:
                    if binid not in bintaxonomy:
                        bintaxonomy[binid] = []
                        bintaxonomy[binid] += [gi2taxid[gi]]
                    else:
                        bintaxonomy[binid] += [gi2taxid[gi]]

        outfile = 'nc_genomes.eukaryotic.tax.txt'
        taxlevels = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies', 'no rank']
        with open(outfile,'w') as out:
            out.write('binid' + '\t' + 'nproteins' + '\t' + '\t'.join(taxlevels) +'\n')
            for binid in bintaxonomy:
                taxids = bintaxonomy[binid]
                if len(taxids) >= 2:
                    taxobjects = [ taxopy.Taxon(taxid,taxdb) for taxid in taxids] 
                    major_vote = taxopy.find_majority_vote(taxobjects, taxdb)
                    major_taxid = major_vote.taxid
                    lineage = lineage_to_name(major_taxid,taxdb)
                    taxstring = '\t'.join([lineage[l] for l in taxlevels])
                    nproteins = len(taxids)
                    out.write(binid + '\t' + str(nproteins) + '\t' + taxstring + '\n')

                else:
                    lineage = lineage_to_name(taxids[0],taxdb)
                    taxstring = '\t'.join([lineage[l] for l in taxlevels])
                    nproteins = 1
                    out.write(binid + '\t' + str(nproteins) + '\t' + taxstring + '\n')






    def assign_uniprot_taxonomy(args,taxdb,bins,vogclades):    

        
        vogclades_tax = {k:None for k in vogclades}
        for cladeid in vogclades:
            if cladeid in name2taxid:
                taxid = name2taxid[cladeid]
                lineage = lineage_to_name(taxid,taxdb)
                lineagename = [ lineage[lin] for lin in ['superkingdom','phylum','class','order','family','genus','species','subspecies'] ]
                vogclades_tax[cladeid] = lineagename
            elif cladeid == 'CressDNAParvo':
                lineagename = ['Viruses','NA','NA','CRESS','CRESS','NA','NA','NA']
                vogclades_tax[cladeid] = lineagename
            elif cladeid == 'PolyoPapillo':
                taxid = name2taxid['Papillomaviridae']
                lineage = lineage_to_name(taxid,taxdb)
                lineagename = [ lineage[lin] for lin in ['superkingdom','phylum','class','order','family','genus','species','subspecies'] ]
                vogclades_tax[cladeid] = lineagename
            elif cladeid == 'NCLDV':
                lineagename = ['Viruses','NA','NA','NCLDV','NCLDV','NA','NA','NA']
                vogclades_tax[cladeid] = lineagename


        ### Read bins that are clearly crass-like
        crassbins = read_crassphage_hits(args)

        ### Uniprot
        uniprot_tax = parse_uniprot_hits(args,taxdb,bins)

        ### rvdb / should be primarily eukaryotic
        rvdb_tax  = parse_rvdb_hits(args,taxdb,bins)
        
        ### Summarise taxonomy for each binid 
        crasstaxonomy = ['Viruses','NA','NA','Caudovirales','crAss-like','unclassified Podoviridae','crAss-like viruses','NA']

        fileout = os.path.join(args.v, 'checkv_refined_taxonomy.txt')
        with open(fileout,'w') as out:
            out.write('binid\tcompleteness\trefined_taxonomy\tClusterTaxonomy\tVOGclade\n')
            for binid in bins:
                binobject = bins[binid]
                unitax = uniprot_tax[binid]
                rvdbtax = rvdb_tax[binid]
                vog_cluster_taxonomy = bins[binid].clustertax
                VOGclade = bins[binid].VOGclade

                ### We're confident in Crass bins
                if binid in crassbins:
                    refined_tax = ';'.join(crasstaxonomy)
                    VOGclade = 'CrAss-Like'
                    lineout = [binid,binobject.completeness,refined_tax,refined_tax,VOGclade]
                    lineout = '\t'.join(lineout) + '\n'
                    out.write(lineout)
                else:

                    ### Compare the VOGclade with the taxonomy of RVDB, Uniprot and VOG
                    refined_tax = None
                    for tax in [rvdbtax,unitax,vog_cluster_taxonomy]:
                        if VOGclade in tax:
                            refined_tax = tax
                    if refined_tax is None and VOGclade != 'NA':
                        refined_tax = ';'.join(vogclades_tax[VOGclade])
                    elif refined_tax is None:
                        refined_tax = vog_cluster_taxonomy
                    
                    lineout = [binid,binobject.completeness,refined_tax,vog_cluster_taxonomy,VOGclade]
                    lineout = '\t'.join(lineout) + '\n'
                    out.write(lineout)



    ### Read bins and their taxonomy based on the VOG
    bins,vogclades = read_checkv_quality(args)

    ### 
    assign_uniprot_taxonomy(args,taxdb,bins,vogclades)

   
### Orthologous proteins for making Phylogenetic Tree
def ortholog_proteins(args):
    



        def get_terl_proteins(checkv_directory):
            '''
            VOG terminase markers are found with a simple grep-command: " grep -i terminase vog.annotations.tsv  | grep -i large  " 
            '''
            terl_proteins = set()
            terl_genomes = dict()
            
            ### Annotation ID's 

            VOG = set(['VOG00419','VOG00699','VOG00709','VOG00731','VOG00732','VOG01032','VOG01094','VOG01180','VOG01426'])
            crassterminase = 'YP_009052554.1'


            ### Annotation files to parse 
            blastp_crass = os.path.join(checkv_directory,'nc_genomes.crassphage.polterm.m6')
            voghmmfile = os.path.join(checkv_directory,'nc_genomes.VOG.hmmsearch')
            eggfile = os.path.join(checkv_directory,'eggnog/nc_genomes_proteins.emapper.annotations')
            interprofile = os.path.join(checkv_directory,'eggnog/interproscan.tsv')

            ### Crassphage markers
            with open(blastp_crass,'r') as infile:
                for line in infile:
                    line = line.strip().split('\t')
                    proteinid = line[0]
                    target = line[1]
                    binid = '_'.join(proteinid.split('_')[:2])
                    aln = int(line[3])
                    evalue = float(line[10])
                    if aln < 350:
                        continue
                    if target == crassterminase:
                        terl_proteins.add(proteinid)
                        if not binid in terl_genomes:
                            terl_genomes[binid] = set([target])
                        else:
                            terl_genomes[binid].add(target)
            
            ### VOG markers 
            with open(voghmmfile,'r') as infile:
                for line in infile:
                    if line[0] =='#':
                        continue
                    line = line.strip().split(' ')
                    line = [ x for x in line if x != '' ]
                    proteinid = line[0]
                    binid = '_'.join(proteinid.split('_')[:2])
                    vog = line[2]
                    score = float(line[5])
                    if score >= 30 and vog in VOG:
                        terl_proteins.add(proteinid)
                        if not binid in terl_genomes:
                            terl_genomes[binid] = set([vog])
                        else:
                            terl_genomes[binid].add(vog)

            ### Fluffy String Search for Interpro and Eggnog
            ### Eggnog
            with open(eggfile,'r') as infile:
                for line in infile:
                    if line[0] == '#':
                        continue
                    line = line.strip().split('\t')
                    proteinid = line[0]
                    egg = line[1]
                    binid = '_'.join(proteinid.split('_')[:2])
                    eggdescription = line[-1].lower()
                    if 'terminase' in eggdescription and 'large' in eggdescription:
                        terl_proteins.add(proteinid)
                        if not binid in terl_genomes:
                            terl_genomes[binid] = set([egg])
                        else:
                            terl_genomes[binid].add(egg)
            ### Interpro 
            with open(interprofile,'r') as infile:
                for line in infile:
                    line = line.strip().split('\t')
                    proteinid = line[0] 
                    binid = '_'.join(proteinid.split('_')[:2])
                    IPRid = line[-2]
                    iprdesc = line[-1].lower()
                    if 'terminase' in iprdesc and 'large' in iprdesc:
                        terl_proteins.add(proteinid)
                        if not binid in terl_genomes:
                            terl_genomes[binid] = set([IPRid])
                        else:
                            terl_genomes[binid].add(IPRid)
            return terl_proteins, terl_genomes

        def write_out_proteins(checkv_directory,marker_proteins,genome_taxonomy,marker,viralfamily=None):
            print(viralfamily)
            prodigalfile_AA = os.path.join(checkv_directory,'nc_genomes_proteins.faa')
            if not viralfamily is None:
                protein_outfile = os.path.join(checkv_directory,'orthologs',marker +'.'+ viralfamily + '.faa')
                written_proteins = set()
                outhandle = open(protein_outfile , 'w')
                for record in SeqIO.parse(open(prodigalfile_AA, 'r'), 'fasta'):
                    recordname = record.id
                    genomeid = '_'.join(recordname.split('_')[:2])
                    if recordname in marker_proteins and genomeid in genome_taxonomy:
                        if recordname in written_proteins:
                            continue
                        recordtaxonomy = genome_taxonomy[genomeid][0]
                        if recordtaxonomy == viralfamily:
                            written_proteins.add(recordname)
                            record.description = ''
                            SeqIO.write(record, outhandle, 'fasta')
                
            else:
                protein_outfile = os.path.join(checkv_directory,'orthologs',marker + '.faa')
                written_proteins = set()
                i = 0
                print(protein_outfile)
                outhandle = open(protein_outfile , 'w')
                for record in SeqIO.parse(open(prodigalfile_AA, 'r'), 'fasta'):
                    recordname = record.id
                    if recordname in marker_proteins:
                        if recordname in written_proteins:
                            continue
                        written_proteins.add(recordname)
                        SeqIO.write(record, outhandle, 'fasta')
                        i += 1

        def make_phylo_Tree(checkv_directory,marker,viralfamily=None):

            '''
            load the following: mafft/7.453 trimal/1.4.1 iqtree/1.6.8
             
            Prepare Phylofile 
            Run IQtree
            '''
            if not viralfamily is None:
                protein_outfile = os.path.join(checkv_directory,'orthologs',marker + '.'  + viralfamily +'.faa')
            else:
                protein_outfile = os.path.join(checkv_directory,'orthologs',marker + '.faa')

            alnfile = protein_outfile.replace('faa','aln')
            
            ### Run Mafft
            try:
                mafftexecutable = '/services/tools/mafft/7.453/bin/mafft'
                with open(alnfile,'w') as outfile:
                    command = [mafftexecutable,
                            '--auto', protein_outfile]
                    subprocess.run(command, stdout= outfile)
            except:
                message = 'ERROR: MAFFT finished abnormally.'
                print(message)
                sys.exit(1)


            ### Exchange X to - 
            # This command screws with the protein names... 
            #alnfilex = alnfile.replace('aln','xaln')
            #ps = subprocess.Popen(['cat',alnfile], stdout = subprocess.PIPE)
            #with open(alnfilex,'w') as outfile:
            #    command = ['tr',''''X'''',"'-'"]
            #    subprocess.run( command, stdin=ps.stdout, stdout=outfile)

            ### Run Trimal 
            try:
                trimalexecutable = '/services/tools/trimal/1.4.1/bin/trimal'
                phyfile = protein_outfile.replace('faa','phy')
                command = [trimalexecutable,
                        '-in', alnfile,
                        '-phylip',
                        '-out', phyfile]
                subprocess.check_call(command)
            except:
                message = 'ERROR: Trimal finished abnormally.'
                print(message)
                sys.exit(1)

            # module load iqtree/1.6.8
            # iqtree -s terl.crAss-like.phy -m TEST -nt 24 -safe -bb 1000 -alrt 1000 -redo
            # Best-fit model: VT+F+G4 chosen according to BIC
            # iqtree -s terl.phy -m LG -nt 24
            #try:
            #    iqtree_executable = '/services/tools/iqtree/1.6.8/bin/iqtree'
            #    command = [iqtree_executable,
            #            '-iqtree',
            #            '-s', phyfile,
            #            '-m', 'LG',
            #            '-nt','24']
            #    subprocess.check_call(command)
            #except:
            #    message = 'ERROR: IQtree finished abnormally.'
            #    print(message)
            #    sys.exit(1)

        def write_tree_annotation(checkv_directory,marker,terl_genomes,terl_proteins):

            tree_annotation_file = os.path.join(checkv_directory,'orthologs',marker + '.annotation.txt')


            
            with open(tree_annotation_file,'w') as out:
                out.write('proteinid\tgenomeid\tcluster\tfamily\tVOGclade\n')
                for proteinid in terl_proteins:
                    genomeid = '_'.join(proteinid.split('_')[:2])
                    family, VOGclade = genome_taxonomy[genomeid]
                    cluster = genomeid.split('_')[1]
                    out.write('{}\t{}\t{}\t{}\t{}\n'.format(proteinid,genomeid,cluster,family,VOGclade))
            
        
        

        
        ### Get all Putative Terminase (large subunit) proteins 
        checkv_directory = args.v
        terl_proteins, terl_genomes = get_terl_proteins(checkv_directory)
        
        ### Load Taxonomy
        genome_taxonomy = dict()
        with open(os.path.join(checkv_directory,'checkv_refined_taxonomy.txt') ,'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                genomeid = line[0]
                if genomeid in terl_genomes:
                    tax = line[2]
                    VOGclade = line[4]
                    family = tax.split(';')[4]
                    genome_taxonomy[genomeid] = [family,VOGclade]


        write_out_proteins(checkv_directory,terl_proteins,genome_taxonomy,'terl',viralfamily='crAss-like')
        make_phylo_Tree(checkv_directory, marker = 'terl',viralfamily='crAss-like')
        write_tree_annotation(checkv_directory,  'terl', terl_genomes, terl_proteins)




if __name__ == "__main__":
    args = parser.parse_args()

    ### Write filtered quality summary file
    checkv_parsers.write_filtered_quality_summary(args)
    
    nc_viral_bins_file = os.path.join(args.v,'nc_viralbins.txt')
    if os.path.exists(nc_viral_bins_file):
        print('Do you really wanna rerun the whole thing again? Remove this file to make it again:',nc_viral_bins_file)
    else:
        ### Parse checkV files and calculate VOG taxonomy for each Bin
        parse_checkv_write_overview(args)

    ####### Extensive Annotation of NC bins #######
    ncbins = set()
    with open(nc_viral_bins_file,'r') as infile:
        for line in infile:
            binid = line.strip()
            ncbins.add(binid)

    ### Subset the Cleaned Genome fasta with the >80% complete bins
    # Annotate Proteins to do more Extensive Taxonomy Characterisation
    if len(ncbins) > 0:

        checkv_directory = args.v
        prodigalfile_AA = os.path.join(checkv_directory,'nc_genomes_proteins.faa')
        if os.path.exists(prodigalfile_AA):
            print('Proteins of NC genomes already predicted, delete first to recreate ', prodigalfile_AA)
        else:
            tax_annotations.write_out_proteomes(args,ncbins,prodigalfile_AA)

        ###  Annotation of Proteomes
        annotation_ncgenomes_proteomes(args)

        fileout = os.path.join(args.v, 'checkv_refined_taxonomy.txt')
        if os.path.exists(fileout):
            print('The ',fileout, 'already exists. Delete it to recreate it!')
        else:
            write_extended_taxonomy(args)

    else:
        print('No NC genomes found ???')
        sys.exit(1)


