#!/bin/python

import os 
import sys
from Bio import SeqIO


def load_checkv_db_files(checkv_db_dir):
    
    ### Load GCA dir.
    gcatax = {}
    checkv_gca_taxonomy = os.path.join(checkv_db_dir,'S2_taxonomy.csv')
    with open(checkv_gca_taxonomy,'r') as infile:
        header = infile.readline().strip().split(';')
        lineage_index = header.index('NCBI lineage')
        gca_size_index = header.index('Genome length')
        vog_clade_index = header.index('VOG clade')
        for line in infile:
            line = line.strip().split(';') 
            gca = line[0]
            if gca not in gcatax:
                gcatax[gca] = {'size':None,'vogclade':None,'lineage':None}
                size = line[gca_size_index]
                vogclade = line[vog_clade_index]
                lineage = line[lineage_index:]
                lineage = [x.replace('"','').replace(' ','') for x in lineage]
                lineage = ','.join(lineage)
                gcatax[gca]['size'] = size
                gcatax[gca]['vogclade'] = vogclade
                gcatax[gca]['lineage'] = lineage
    
    ### Rep tax 
    DTRtax = {}
    vogclades = set()
    checkv_dtr_taxonomy = os.path.join(checkv_db_dir,'S3_refs.tsv')
    with open(checkv_dtr_taxonomy,'r') as infile:
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
    checkv_clusters = os.path.join(checkv_db_dir,'S4_cluster.csv')
    with open(checkv_clusters,'r') as infile:
        header = infile.readline().strip().split(';')
        gca_index = header.index('Rep. GenBank')
        for line in infile:
            line = line.strip().split(';')
            gca_rep = line[gca_index]
            if gca_rep == '-':
                gca_rep = None
            genomes = line[1].split(',')
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

            binobject = bins[binid]
            quality = line[quality_index]
            genome_copies = float(line[genome_copies_index])
            contaminaton = line[contamination_index]
            
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
                    if aai_confidence == 'high':
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

    def transform_lineage(lineage,VOGclade):
        '''
        To make a consistency between the Taxonomy of checkV lineage strings and VOG taxonomy,
                make a Tax-string of 8 nodes for Caudovirales

        '''

        taxlevels = [[0,'superkingdom'],
                [1,'phylum'],
                [2,'class'],
                [3,'order'],
                [4,'family'],
                [5,'genus'],
                [6,'species'],
                [7,'strain']]

        if VOGclade == 'Caudovirales':
            lineage_split = lineage.split(',')
            new_lineage = [ lineage_split[0] , 'NA', 'NA' ] 
            new_lineage += lineage_split[1:]
            fill_NAs = [ 'NA' for i in range(8-len(new_lineage))]
            new_lineage += fill_NAs
        elif VOGclade == 'NCLDV':
            lineage_split = lineage.split(',')
            new_lineage = [ lineage_split[0] , 'NA', 'NA' ] 
            new_lineage += lineage_split[1:]
            fill_NAs = [ 'NA' for i in range(8-len(new_lineage))]
            new_lineage += fill_NAs
        else:
            lineage_split = lineage.split(',')
            new_lineage = lineage.split(',')
            new_lineage = [ lineage_split[0] , 'NA', 'NA' ] 
            new_lineage += lineage_split[1:]
            fill_NAs = [ 'NA' for i in range(8-len(new_lineage))]
            new_lineage += fill_NAs

        return(new_lineage)



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
                lineage = motherbins[motherbin].motherbintax_GCA
                VOGclade =  motherbins[motherbin].motherbintax_VOGclade
                clustertaxonomy = transform_lineage(lineage,VOGclade) 
                clustertaxonomy = ';'.join(clustertaxonomy)

            if not binobject.GCAtax is None:
                lineage = binobject.GCAtax 
                Taxonomy = transform_lineage(lineage,VOGclade) 
                Taxonomy = ';'.join(Taxonomy)
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
        for item in ncbins:
            out.write(str(item)+'\n')


    ### Parse checkV generated fasta with Host-contamination removed 
    checkv_directory = args.v

    fasta_file  = os.path.join(checkv_directory,'cleaned_contigs.fna')
    ncgenomes_outfile = os.path.join(checkv_directory,'nc_genomes.fna')
    written_bin_records = set()
    outhandle = open(ncgenomes_outfile , 'w')
    for record in SeqIO.parse(open(fasta_file, 'r'), 'fasta'):
        recordname = record.id 
        binid = '_'.join(recordname.split('_')[0:2])
        if binid in ncbins:
            if not binid in written_bin_records:
                new_record = binid
                record.id = new_record 
                record.description = ''
                SeqIO.write(record, outhandle, 'fasta')
                written_bin_records.add(binid)
