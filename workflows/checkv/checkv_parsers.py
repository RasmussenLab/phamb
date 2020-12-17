#!/bin/python

import os 
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import taxopy



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

def load_checkv_db_files(checkv_db_dir,name2taxid,taxdb):
    '''
    Load and organise Taxonomy of all Reference Genomes in the checkV database. 
    '''
    
    lineage_levels = ['superkingdom','phylum','class','order','family','genus','species','subspecies']

    ### Load GCA dir.
    gcatax = {}
    checkv_gca_taxonomy = os.path.join(checkv_db_dir,'checkv_genbank.tsv')
    with open(checkv_gca_taxonomy,'r') as infile:
        header = infile.readline().strip().split('\t')
        lineage_index = header.index('lineage')
        gca_size_index = header.index('genome_length')
        vog_clade_index = header.index('vog_clade')
        for line in infile:
            line = line.strip().split('\t') 
            gca = line[0]
            if gca not in gcatax:
                gcatax[gca] = {'size':None,'vogclade':None,'lineage':None}
                size = line[gca_size_index]
                vogclade = line[vog_clade_index]
                lineage = line[lineage_index].replace(' ','').split(';')

                ### Make complete Lineage
                lowest_taxid = None
                while lowest_taxid == None:
                    name = lineage.pop(-1)
                    if name in name2taxid:
                        lowest_taxid = name2taxid[name]                
                lineage = lineage_to_name(lowest_taxid,taxdb) 
                taxstring = ';'.join([lineage[l] for l in lineage_levels])
                gcatax[gca]['size'] = size
                gcatax[gca]['vogclade'] = vogclade
                gcatax[gca]['lineage'] = taxstring
    
    ### DTR taxonomy 
    DTRtax = {}
    vogclades = set()
    checkv_DTR_taxonomy = os.path.join(checkv_db_dir,'checkv_circular.tsv')
    with open(checkv_DTR_taxonomy,'r') as infile:
        header = infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            checkvid = line[0]
            vogclade = line[5]
            vogclades.add(vogclade)
            if checkvid not in DTRtax:
                DTRtax[checkvid] = vogclade


    ### establish Cluster network with Representative Genomes 
    checkv_genome_gcatax = {}
    checkv_clusters = os.path.join(checkv_db_dir,'checkv_clusters.tsv')
    with open(checkv_clusters,'r') as infile:
        header = infile.readline().strip().split('\t')
        gca_index = header.index('genbank_rep')
        cluster_members = header.index('genome_ids')
        for line in infile:
            line = line.strip().split('\t')
            gca_rep = line[gca_index]
            if gca_rep == 'NULL':
                gca_rep = None
            genomes = line[cluster_members].split(',')
            for genome in genomes:
                checkv_genome_gcatax[genome] = gca_rep

    return checkv_genome_gcatax, gcatax, DTRtax


def read_checkv_quality(args,bins):
    '''
    Load annotation of each Bin from the checkV quality_summary file 
    '''

    quality_file = os.path.join(args.v,'quality_summary.tsv')

    with open(quality_file,'r') as infile:
        header = infile.readline().strip().split('\t')
        genome_copies_index = header.index('genome_copies')
        #contamination_index = header.index('contamination')
        quality_index = header.index('checkv_quality')
        method_index = header.index('completeness_method')
        miuvig_quality = header.index('miuvig_quality')
        viral_genes_index = header.index('viral_genes')
        prophage_index = header.index('prophage')
        completeness_index = header.index('completeness')

        for line in infile:
            line = line.strip().split('\t')
            binid = line[0]

            binobject = bins[binid]
            quality = line[quality_index]
            genome_copies = float(line[genome_copies_index])
            #contaminaton = line[contamination_index]
            
            ### Filter off Bins with little Viral evidence
            if quality == 'Not-determined' and line[viral_genes_index] == '0':
                continue
            if genome_copies > 1.25:
                continue
            
            binobject.checkv_results['completeness_method'] = line[method_index]
            binobject.checkv_results['completeness'] = line[completeness_index]
            binobject.checkv_results['prophage'] = line[prophage_index]
            binobject.checkv_results['checkv_quality'] = quality
            binobject.checkv_results['miuvig_quality'] = line[miuvig_quality]
            bins[binid] = binobject
    return bins


def label_checkv_bins_by_quality(args):
    """[summary]

    Args:
        args ([type]): [description]
        bins ([type]): [description]
    """
    quality_file = os.path.join(args.v,'quality_summary.tsv')

    population_type_score = {3:'HQ-ref',2:'Grey-matter',1:'Dark-matter'}

    viral_population_types = defaultdict()
    with open(quality_file,'r') as infile:
        header = infile.readline().strip().split('\t')
        genome_copies_index = header.index('genome_copies')
        contamination_index = header.index('contamination')
        quality_index = header.index('checkv_quality')
        method_index = header.index('completeness_method')
        miuvig_quality_index = header.index('miuvig_quality')
        viral_genes_index = header.index('viral_genes')
        prophage_index = header.index('prophage')
        completeness_index = header.index('completeness')

        for line in infile:
            line = line.strip().split('\t')
            binid = line[0]
            motherbin = binid.split('_')[1]
            population_type = None

            genome_copies = float(line[genome_copies_index])
            contaminaton = line[contamination_index]
            
            completeness_method = line[method_index]
            completeness = line[completeness_index] 
            checkv_quality = line[quality_index]
            miuvig_quality = line[miuvig_quality_index]

            if genome_copies > 1.25:
                continue

            if checkv_quality in ['High-quality'] and completeness_method == 'AAI-based':
                population_type = 3
                if motherbin in viral_population_types:
                    if population_type > viral_population_types[motherbin]:
                        viral_population_types[motherbin] = population_type
                else:
                    viral_population_types[motherbin] = population_type
                continue
            
            if checkv_quality in ['High-quality','Medium-quality'] and completeness_method == 'HMM-based':
                population_type = 2
                if motherbin in viral_population_types:
                    if population_type > viral_population_types[motherbin]:
                        viral_population_types[motherbin] = population_type
                else:
                    viral_population_types[motherbin] = population_type
                continue
            
            if checkv_quality in ['Medium-quality']:
                population_type = 2
                if motherbin in viral_population_types:
                    if population_type > viral_population_types[motherbin]:
                        viral_population_types[motherbin] = population_type
                else:
                    viral_population_types[motherbin] = population_type
                continue

            
            population_type = 1
            if not motherbin in viral_population_types:
                viral_population_types[motherbin] = population_type
        
        fileout = os.path.join(args.v,'population_types.tsv')
        with open(fileout,'w') as out:
            for motherbin in viral_population_types:
                _type = viral_population_types[motherbin]
                population_type = population_type_score[_type]
                lineout = '\t'.join( [motherbin,population_type] )
                out.write(lineout + '\n')




    



def write_filtered_quality_summary(args):
    '''
    Write a filtered edition of the quality file for those genomes without a viral region (only viral genes between host genes)
    '''

    contamination_file = os.path.join(args.v,'contamination.tsv')
    safe_viral_ids = set() 
    with open(contamination_file,'r') as infile:
        header = infile.readline().strip().split('\t')
        viral_length_index = header.index('viral_length')
        viral_gene_index = header.index('viral_genes')
        region_index = header.index('region_types')
        for line in infile:
            line = line.strip().split('\t')
            genomeid = line[0]
            regiontype = line[region_index]
            if int(line[viral_length_index]) != 0 or int(line[viral_gene_index]) != 0:
                safe_viral_ids.add(genomeid)
            if regiontype == 'unclassified':
                safe_viral_ids.add(genomeid)
    
    quality_summary_file = os.path.join(args.v,'quality_summary.tsv')
    quality_summary_file_filtered = os.path.join(args.v,'quality_summary_filtered.tsv')

    with open(quality_summary_file,'r') as infile, open(quality_summary_file_filtered,'w') as outfile:
        header = infile.readline()
        outfile.write(header)
        for line in quality_summary_file:
            genomeid = line.strip().split('\t')[0]
            if genomeid in safe_viral_ids:
                outfile.write(line)
                

                






def read_checkv_completeness(args,bins,motherbins,checkv_genome_gcatax,gcatax,DTRtax):
    '''
    Load annotation of each Bin from the checkV completeness file 
    '''

    completeness_file = os.path.join(args.v,'completeness.tsv')


    with open(completeness_file,'r') as infile:
        header = infile.readline().strip().split('\t')
        aai_confidence_index = header.index('aai_confidence')
        aai_tophit_index = header.index('aai_top_hit')
        hmm_name_index = header.index('hmm_name')
        for line in infile:
            line = line.strip().split('\t')
            binid = line[0]
            motherbinid = binid.split('_')[1]
            binobject = bins[binid]

            ### checkV reference genome used - look for GCA taxonomy
            if binobject.checkv_results['completeness_method'] == 'AAI-based':
                checkv_reference = line[aai_tophit_index]
                aai_confidence = line[aai_confidence_index]
                binobject.checkv_reference = checkv_reference
                binobject.checkv_AAI_confidence = aai_confidence
                if checkv_reference in DTRtax:
                    vogclade = DTRtax[checkv_reference] 
                binobject.VOGclade = vogclade

                ### Lookup 
                gca_rep = checkv_genome_gcatax[checkv_reference]
                if not gca_rep is None:
                    lineage = gcatax[gca_rep]['lineage']
                    binobject.GCAtax = lineage
                    vogclade = gcatax[gca_rep]['vogclade']

                    ### Save for GCA taxonomy for the whole Cluster of bins - if the confidence is high 
                    if aai_confidence in ['medium','high']:
                        clusterobject = motherbins[motherbinid]
                        clusterobject.motherbintax_GCA = lineage
                        clusterobject.motherbintax_VOGclade = vogclade
                        motherbins[motherbinid] = clusterobject

            else:
                hmm_profile_name = line[hmm_name_index]
                binobject.checkv_results['HMMname'] = hmm_profile_name
            
            bins[binid] = binobject
    return bins, motherbins



def write_bin_overview(args,bins,motherbins):
    '''
    Write out some Taxonomy info for each Bin 
    - If HMM classification has been used or no GCA ref available - use VOG taxonomy 
    - else use Taxonomy of GCA   

    Based on the checkV quality scores for each Bin - selct the NC Viral Genomes and write them to a file
    '''


    ncbins = set()
    outfile = os.path.join(args.v,'checkv_taxonomy.txt')
    with open(outfile,'w') as out:
        out.write('binid\tmotherbin\tcheckv_quality\tcompleteness\tcompleteness_method\tTaxonomy\tClusterTaxonomy\tVOGclade\n')
        for binid in bins:
            binobject = bins[binid]
            VOGclade = binobject.VOGclade
            if binobject.checkv_results['miuvig_quality'] is None:
                continue
            motherbin = binid.split('_')[1]
            
            ### Assign Taxonomy
            Taxonomy = None

            ### Get Taxonomy derived on VOG level
            clustertaxonomy = motherbins[motherbin].motherbintax + ';NA'

            ### If the cluster can be described with a GCA reference with high Confidence use it for all bins ClusterTaxonomy instead of VOG
            if not motherbins[motherbin].motherbintax_GCA is None:
                clustertaxonomy = motherbins[motherbin].motherbintax_GCA
                VOGclade =  motherbins[motherbin].motherbintax_VOGclade


            if not binobject.GCAtax is None:
                Taxonomy = binobject.GCAtax 
            else:
                Taxonomy = binobject.bintax
                Taxonomy = Taxonomy + ';NA'

            ### Extract other relevant metrics for outfile 
            checkv_quality = binobject.checkv_results['checkv_quality']
            completeness = binobject.checkv_results['completeness']
            method = binobject.checkv_results['completeness_method']
            if checkv_quality == 'Complete':
                completeness = '100'

            if completeness != 'NA':
                if float(completeness) >= 80:
                    ncbins.add(binid)
            if VOGclade is None:
                VOGclade = 'NA'
            lineout = '\t'.join( [binid, motherbin,checkv_quality, completeness, method, Taxonomy,clustertaxonomy,VOGclade] )
            out.write(lineout+'\n')
    
    ### Write the names of NC Viral Bins to a file
    nc_viral_bins_file = os.path.join(args.v,'nc_viralbins.txt')
    with open(nc_viral_bins_file,'w') as out:
        for binid in ncbins:
            out.write(str(binid)+'\n')


    ### Parse checkV generated fasta with Host-contamination removed 
    # Note that some of the Viral Genomes are divided into fragments due to the presence of Host-contamination
    # binidX_1 
    # binidX_2 
    # binidX_3
    fragmented_viral_genomes = dict()

    checkv_directory = args.v
    fasta_file  = os.path.join(checkv_directory,'cleaned_contigs.fna')
    ncgenomes_outfile = os.path.join(checkv_directory,'nc_genomes.fna')
    written_bin_records = set()
    outhandle = open(ncgenomes_outfile , 'w')
    for record in SeqIO.parse(open(fasta_file, 'r'), 'fasta'):
        recordname = record.id 
        recordname_split = recordname.split('_')
        if len(recordname_split) == 2:
            binid = '_'.join(recordname.split('_')[0:2])
            if binid in ncbins:
                if not binid in written_bin_records:
                    new_record = binid
                    record.id = new_record 
                    record.description = ''
                    SeqIO.write(record, outhandle, 'fasta')
                    written_bin_records.add(binid)
        else:
            ### Host contamination division.
            binid = '_'.join(recordname.split('_')[0:2])
            if binid in ncbins:
                if not binid in fragmented_viral_genomes:
                    fragmented_viral_genomes[binid] = ""
                    fragmented_viral_genomes[binid] += record.seq
                else:
                    fragmented_viral_genomes[binid] += record.seq

    ### Write the combined fragmenteed viral genomes 
    for binid in fragmented_viral_genomes:
        rec = SeqRecord(fragmented_viral_genomes[binid],id=binid, description='')
        SeqIO.write(rec, outhandle, 'fasta')








def parse_checkv_write_overview(args):

      
    def parse_clusterfile(args):
        '''Return list of contigids, size and affiliation'''

        contigs = {}
        clusterfile = os.path.join(args.c,'clusters.tsv')
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

    nodes = '/home/projects/cpr_10006/projects/phamb/databases/TaxonomyDB/NCBI_Taxonomy_new_db_2019-09-20/nodes.dmp'
    names ='/home/projects/cpr_10006/projects/phamb/databases/TaxonomyDB/NCBI_Taxonomy_new_db_2019-09-20/names.dmp'
    taxdb = taxopy.TaxDb(nodes_dmp=nodes, names_dmp=names, keep_files=True)
    name2taxid = {v:k for k,v in taxdb.taxid2name.items()}
    checkv_db_dir = '/home/projects/cpr_10006/projects/phamb/databases/checkv/checkv-db-v0.6/genome_db'
    checkv_genome_gcatax, gcatax, DTRtax  = load_checkv_db_files(checkv_db_dir,name2taxid,taxdb)
    bins = read_checkv_quality(args,bins)
    bins,motherbins = read_checkv_completeness(args,bins,motherbins,checkv_genome_gcatax,gcatax,DTRtax)

    ### Write final overview of Viral Genomes
    write_bin_overview(args,bins,motherbins)