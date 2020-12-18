#!/bin/python

import sys
import os
import argparse
import gzip
import numpy as np
from collections import defaultdict
from collections import Counter

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
   ''')
parser.add_argument('-v', help='VAMB directory ')
parser.add_argument('-g', help='GTDB-TK taxonomy file')
parser.add_argument('-c', help='CheckV directory')



def get_MAG_overview(args):
    """[Read CheckM file of MQNC MAGs]

    Args:
        args ([type]): [description]
    """
    MAGs = dict()
    checkm_file = '07_binannotation/bacteria/MQNC_MAGS.txt'
    with open(checkm_file,'r') as infile:
        for line in infile:
            line = line.strip()
            genomeid = os.path.basename(line)
            genomeid = genomeid.replace('.fna','')
            MAGs[genomeid] = dict()

    return MAGs

def load_MQ_bins(args):
    ### We cannot include Bins that are far from Complete
    medium_quality_bins = set()
    quality_file = os.path.join(args.c,'quality_summary.tsv')

    if not os.path.exists(quality_file):
        print('Have you run CheckV?')
        sys.exit(1)

    with open(quality_file,'r') as infile:
        header = infile.readline().strip().split('\t')
        genome_copies_index = header.index('genome_copies')
        quality_index = header.index('checkv_quality')
        completeness_index = header.index('completeness')
        
        for line in infile:
            line = line.strip().split('\t')
            binid = line[0]
            motherbin = binid.split('_')[1]
            genome_copies = float(line[genome_copies_index])
            checkv_quality = line[quality_index]

            if genome_copies > 1.25:
                continue

            completeness = line[completeness_index]
            if completeness == 'NA' and checkv_quality == 'Complete':
                completeness = 100
            if completeness == 'NA' and checkv_quality != 'Complete':
                completeness = 0
            if float(completeness) >= 50:
                medium_quality_bins.add(binid)
    return  medium_quality_bins


def load_MAG_taxonomy():
    '''
    GTDB-TK annotation of NC MAGs.
    IDEALLY this should be a consensus vote, not the first occuring...
    If it's annotated all the way to Species level don't parse anymore.
    '''

    gtdb_file = args.g
    lineage = ['superkingdom','phylum','class','order','family','genus','species']
    #dlineage = {lin:None for lin in lineage}

    parsed_motherbins = set()
    MAG_tax = dict()
    with open(gtdb_file,'r') as infile:
        header = infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            motherbin = line[0].split('_')[1]

            if motherbin in MAG_tax:
                if MAG_tax[motherbin]['species'] != 'NA':
                    continue

            MAG_tax[motherbin] = {lin:None for lin in lineage}
            tax = line[2]
            tax = tax.split(';')
            tax = [ x.split('__')[1] for x in tax ]
            for i,lin in enumerate(lineage):
                entry = tax[i]
                if entry == '':
                    entry = 'NA'
                if '_' in entry:
                    entry = entry.split('_')[0]
                if entry == 'Bacteroidota':
                    entry = 'Bacteroidetes'
                
                MAG_tax[motherbin][lin] = entry
    return MAG_tax


def evaluate_host_consistency(_MAG_benchmark):
    """
    1. Load MAG Taxonomy
    2. For each Virus load all Hosts by Spacers 
    3. For each Virus load all Hosts by Prophage hit
    4. For both Methods calculate the most common host at each Lineage 
    5. For both Methods calculate Host-assignment overall Purity
    """ 
    lineage = ['superkingdom','phylum','class','order','family','genus','species']

    ### Load in MAG taxonomy 
    MAG_tax = load_MAG_taxonomy()
    
    ### Calculate the most common Viral host across lineages
    mag_annotations = _MAG_benchmark

    class Virus_class:
        def __init__(self):
            self.domintant_spacer_lineage = None
            self.domintant_prophage_lineage = None
            self.purity_spacer_lineage = None
            self.purity_prophage_lineage = None


    virus_annotations = dict()
    for MAG in mag_annotations:
        if len(mag_annotations[MAG]) == 0:
            continue

        MAG_motherbin = MAG.split('_')[1]

        if not MAG_motherbin in MAG_tax:
            continue
        MAG_lineage = MAG_tax[MAG_motherbin]

        for phage in mag_annotations[MAG]:

            ### 
            if not phage in virus_annotations:
                viral_class = Virus_class()
                viral_class.spacers = {k:[] for k in lineage}
                viral_class.prophages =  {k:[] for k in lineage}
            else:
                viral_class = virus_annotations[phage] 

            host_annotations = mag_annotations[MAG][phage]

            if not host_annotations['spacer'] is None:
                for k in lineage:
                    lin = MAG_lineage[k]
                    viral_class.spacers[k].append(lin) 
            if not host_annotations['prophage'] is None:
                for k in lineage:
                    lin = MAG_lineage[k]
                    viral_class.prophages[k].append(lin) 
            
            virus_annotations[phage] = viral_class
    
    ### Calculate most Common Host Taxonomy at each Lineage 
    ### & Evaluate Host-prediction-purity 
    for phage in virus_annotations:
        viral_class = virus_annotations[phage]

        viral_class.domintant_spacer_lineage = None
        viral_class.domintant_prophage_lineage = None

        dominant_spacer_lineage = {k:'NA' for k in lineage}
        dominant_prophage_lineage = {k:'NA' for k in lineage}
        purity_spacer_lineage = {k:'NA' for k in lineage}
        purity_prophage_lineage = {k:'NA' for k in lineage}

        for k in lineage:
            c = Counter(viral_class.spacers[k])
            cmc = c.most_common()
            if len(cmc) != 0:
                dominant_spacer_lineage[k] = cmc[0][0]
                purity_spacer_lineage[k] = viral_class.spacers[k].count(cmc[0][0])/len(viral_class.spacers[k])

            c = Counter(viral_class.prophages[k])
            cmc = c.most_common()
            if len(cmc) != 0:
                dominant_prophage_lineage[k] = cmc[0][0]
                purity_prophage_lineage[k] = viral_class.prophages[k].count(cmc[0][0])/len(viral_class.prophages[k])
        
        if dominant_spacer_lineage['superkingdom'] == 'NA':
            viral_class.domintant_spacer_lineage = None
        else:
            viral_class.domintant_spacer_lineage = dominant_spacer_lineage
            viral_class.purity_spacer_lineage = purity_spacer_lineage

        if dominant_prophage_lineage['superkingdom'] == 'NA':
            viral_class.domintant_prophage_lineage = None
        else:
            viral_class.domintant_prophage_lineage = dominant_prophage_lineage
            viral_class.purity_prophage_lineage = purity_prophage_lineage
        virus_annotations[phage] = viral_class

    ### Calculate overall Consensus at each Lineage-level for both Methods and calculate consensus scores
    nviruses_both_predictions = {k:None for k in lineage}
    consensus = {k:None for k in lineage}
    for k in lineage:
        results = []
        for phage in virus_annotations:
            viral_class = virus_annotations[phage]
            if not viral_class.domintant_spacer_lineage is None and not viral_class.domintant_prophage_lineage is None:

                if viral_class.domintant_spacer_lineage[k] == viral_class.domintant_prophage_lineage[k]:
                    results.append('consensus')
                else:
                    results.append('not_consesus')
        agreement_percentage = round((results.count('consensus')/len(results))*100,2)
        nviruses_both_predictions[k] = len(results)
        consensus[k] = agreement_percentage
    
    ### Calculate Average Purity 
    overall_spacer_purity = {k:None for k in lineage}
    overall_prophage_purity = {k:None for k in lineage}
    prophage_count = 0 
    spacer_count = 0
    for k in lineage:
        spacer_results = []
        prophage_results = []
        for phage in virus_annotations:
            viral_class = virus_annotations[phage]
            if not viral_class.purity_prophage_lineage is None:
                prophage_results.append( viral_class.purity_prophage_lineage[k] )
            if not viral_class.purity_spacer_lineage is None:
                spacer_results.append( viral_class.purity_spacer_lineage[k] )
        overall_spacer_purity[k] = round(np.mean(spacer_results)*100,4)
        overall_prophage_purity[k] = round(np.mean(prophage_results)*100,4)

        if k == 'superkingdom':
            prophage_count += len(prophage_results)
            spacer_count += len(spacer_results)

    host_taxonomic_results = {'consensus':consensus, 
    'overall_spacer_purity': overall_spacer_purity,
    'overall_prophage_purity': overall_prophage_purity,
    'nviruses_spacer_prediction': spacer_count,
    'nviruses_prophage_prediction': prophage_count,
    'nviruses_both_predictions':nviruses_both_predictions['superkingdom'],
    'total_number_of_viruses':len(virus_annotations)}

    return host_taxonomic_results


def return_host_consistency_table(host_taxonomy_predictions,prediction_method):

    lineage = ['superkingdom','phylum','class','order','family','genus','species']

    header = [[prediction_method,'variable']  + lineage + ['virus_count_method','total_viruses']]
    all_tables = []
    all_tables.append(header)
    for k in host_taxonomy_predictions.keys():
        predictions = host_taxonomy_predictions[k]
        table = []

        row = [prediction_method+'.'+str(k),'spacer_prophage_consensus']
        for lin in lineage:
            row += [predictions['consensus'][lin]]
        row += [predictions['nviruses_both_predictions']]
        row += [predictions['total_number_of_viruses']]
        table.append(row)

        row = [prediction_method+'.'+str(k),'spacer_host_purity']
        for lin in lineage:
            row += [predictions['overall_spacer_purity'][lin]]
        row += [predictions['nviruses_spacer_prediction']]
        row += [predictions['total_number_of_viruses']]
        table.append(row)

        row = [prediction_method+'.'+str(k),'prophage_host_purity']
        for lin in lineage:
            row += [predictions['overall_prophage_purity'][lin]]
        row += [predictions['nviruses_prophage_prediction']]
        row += [predictions['total_number_of_viruses']]
        table.append(row)
        all_tables.append(table)
    return all_tables
        

def contig_to_MAG(args,MAGobject):
    
    clusterfile = os.path.join(args.v,'clusters.tsv')
    CONTIG2MAG = dict()
    with open(clusterfile,'r') as infile:
        for line in infile:
            cluster, contig = line.strip().split()
            sample = contig.split('_')[0]
            MAG = sample + '_' + cluster
            if not MAG in MAGobject:
                continue
            CONTIG2MAG[contig] = MAG
    return CONTIG2MAG




def parse_CRISPR_results(args,MAGs,CONTIG2MAG):

    ### Parse blastn 
    spaceromefile = "07_binannotation/bacteria/snakecrisprcasfinder/new.spacers.viralfragment.m6"
    for f in [spaceromefile]:
        with open(f,'r') as infile:
            for line in infile:
                line = line.strip().split()
                spacername = line[0]
                contigname = spacername.split(':')[0]

                ### Watch out for Proteein number Suffix added to each contigname!
                contigname = '_'.join(contigname.split('_')[:-1])
                MAG = CONTIG2MAG[contigname]

                MAGmotherbin = MAG.split('_')[1]
                phagebin = line[1]
                sample = phagebin.split('_')[0]
                phagemotherbin = phagebin.split('_')[1]
                identity = float(line[2])
                spacer_covered,spacer_len = int(line[3]),int(line[-2])
                bitscore = float(line[-3])
                evalue = float(line[-4])
                mismatches = int(line[4])
                
                ### Stringent filtering of hits
                if mismatches > 2:
                    continue
                if identity >= 95 and spacer_covered/spacer_len >=0.95:

                    if not phagemotherbin in MAGs[MAG]:
                        MAGs[MAG][phagemotherbin] = {'spacer':True,'prophage':None,'both':None}
    return MAGs




def parse_FastANI_results(args,MAGs,medium_quality_bins,population_type,number_of_5000bp_pieces):
    
    print('Parsing FastANI hits')
    prophage_file = os.path.join("08_crisprcas",'fastani','MGXVIR.fastani.txt.gz')
    with gzip.open(prophage_file,'rt') as infile:
        for line in infile:
            line = line.strip().split('\t')
            ANI = float(line[2])
            coverage = (int(line[3])/int(line[4]))*100
            n = int(line[3])
            if ANI >= 90 and n >= number_of_5000bp_pieces or coverage >= 75:
                phagebin = os.path.basename(line[0]).replace('.fna','')

                if not phagebin in medium_quality_bins:
                    continue

                sample = phagebin.split('_')[0]
                phagemotherbin = phagebin.split('_')[1]
                if not str(phagemotherbin) in population_type:
                    continue
                ctype = population_type[str(phagemotherbin)]

                MAG = os.path.basename(line[1]).replace('.fna','')
                bacterial_clusters = MAG.split('_')[1]

                if not phagemotherbin in MAGs[MAG]:
                    MAGs[MAG][phagemotherbin] = {'spacer':None,'prophage':True,'both':None,'coverage':coverage}
                else:
                    if phagemotherbin in MAGs[MAG] and MAGs[MAG][phagemotherbin]['spacer'] is True:
                        MAGs[MAG][phagemotherbin] = {'spacer':True,'prophage':True,'both':True,'coverage':coverage}
    return MAGs


def parse_CRISPR_PROPHAGE_hits(args):
    """Function for parsing CRISPR-spacer blast-files

    Args:
        args (arguments from commandline): [arguments from commandline]

    """
    ### Parse MAG spacers
    _MAG_benchmark = {k:get_MAG_overview(args) for k in range(1,2) }
    
    MAGs_by_sample = get_MAG_overview(args)
    CONTIG2MAG = contig_to_MAG(args, _MAG_benchmark[1] )


    ### Read in Viral population types
    population_type = dict()
    with open(os.path.join(args.c,'population_types.tsv'),'r') as infile:
        for line in infile:
            phagemotherbin, ctype = line.strip().split()
            population_type[phagemotherbin] = ctype

    medium_quality_bins = load_MQ_bins(args)
    
    for i in _MAG_benchmark:
        _MAG_benchmark[i] = parse_CRISPR_results(args,_MAG_benchmark[i],CONTIG2MAG)
    _MAG_benchmark_darkmatter = _MAG_benchmark.copy()

    ### Parse Prophage hits for HQ-Viruses
    for i in _MAG_benchmark:
        _MAG_benchmark[i] = parse_FastANI_results(args, _MAG_benchmark[i], medium_quality_bins,population_type, number_of_5000bp_pieces= i)

    ### Calculate Host-prediction accuracy  
    host_results = {i:None for i in _MAG_benchmark}
    for i in _MAG_benchmark:
        host_results[i] =  evaluate_host_consistency(_MAG_benchmark[i]) 

    host_prediction_results = return_host_consistency_table(host_results,prediction_method='fastANI')
    combined_table = np.array(host_prediction_results)
    combined_table = np.concatenate(combined_table,axis=0)

    if not os.path.exists( os.path.join("08_crisprcas",'crispr_prophage_results')):
        os.makedirs(os.path.join('08_crisprcas','crispr_prophage_results'))

    fileout = os.path.join('08_crisprcas','crispr_prophage_results','FASTANI_host_predictions.tsv')
    np.savetxt(fileout, combined_table , delimiter='\t', header='', newline='\n', fmt='%s')


    ### Write results
    results_directory = os.path.join("08_crisprcas",'crispr_prophage_results')

    if not os.path.exists(results_directory):
        os.makedirs(results_directory)

    for i in _MAG_benchmark:
        MAGs = _MAG_benchmark[i]
        bacterial_genome_decision = {MAG:{'spacer':None,'prophage':None,'both':None} for MAG in MAGs}
        outfile = os.path.join(results_directory,'MAG_crispr_prophage_summary.'+str(i*5000)+'.txt')
        with open(outfile,'w') as out:
            for MAG in MAGs:
                bacterial_cluster = MAG.split('_')[1]
                if len(MAGs[MAG]) == 0:
                    lineout = '\t'.join( [bacterial_cluster,MAG,'NA','NA','NA','NA'])
                    out.write(lineout+'\n')
                else:
                    for phagemotherbin in MAGs[MAG]:
                        MAG_sample = MAG.split('_')[0]


                        spacer = MAGs[MAG][phagemotherbin]['spacer']
                        prophage =  MAGs[MAG][phagemotherbin]['prophage']
                        both = MAGs[MAG][phagemotherbin]['both']
                        if 'coverage' in MAGs[MAG][phagemotherbin]:
                            coverage = round(MAGs[MAG][phagemotherbin]['coverage'],2)
                        else:
                            coverage = 0
                        if not str(phagemotherbin) in population_type:
                            continue

                        ctype = population_type[str(phagemotherbin)]

                        if both is True: 
                            viral_association = 'both'
                        else:
                            if spacer is True:
                                viral_association = 'spacer'
                            else:
                                viral_association = 'prophage'

                    
                        lineout = '\t'.join( [bacterial_cluster,MAG,phagemotherbin,ctype,viral_association, str(coverage)] )
                        out.write(lineout+'\n')

                        ### Only make decisions based on HQ-ref Viruses
                        if ctype in ['HQ-ref']:
                            bacterial_genome_decision[MAG][viral_association] = True

        outfile = os.path.join(results_directory,'MAG_crispr_prophage_decision.'+str(i*5000)+'.txt')
        with open(outfile,'w') as out:
            for MAG in bacterial_genome_decision:
                bacterial_cluster = MAG.split('_')[1]
                viral_association = 'None'
                spacer = bacterial_genome_decision[MAG]['spacer']
                prophage =  bacterial_genome_decision[MAG]['prophage']
                both = bacterial_genome_decision[MAG]['both']

                if prophage is True and spacer is True:
                    viral_association = 'both'
                    lineout = '\t'.join( [bacterial_cluster,MAG,viral_association] )
                    out.write(lineout+'\n')
                else:
                    if prophage is True:
                        viral_association = 'prophage'
                    elif spacer is True:
                        viral_association = 'spacer'
                    lineout = '\t'.join( [bacterial_cluster,MAG,viral_association] )
                    out.write(lineout+'\n')




def parse_CRISPR_PROPHAGE_hits_blastn(args):
    '''
    Annotate integrated/aliggned phages similar to Nature Paper: https://www.nature.com/articles/s41587-020-0718-6#Sec9 
    '''

    phagesize_ratios = [1,1.5,2]
    _MAG_benchmark_blastn = {k:get_MAG_overview(args) for k in range(len(phagesize_ratios)) }
    

    ### Read in Viral population types
    population_type = dict()
    with open(os.path.join(args.c,'population_types.tsv'),'r') as infile:
        for line in infile:
            phagemotherbin, ctype = line.strip().split()
            population_type[phagemotherbin] = ctype
    
    ### Parse cluster-file
    CONTIG2MAG = contig_to_MAG(args, _MAG_benchmark_blastn[0] )
    medium_quality_bins = load_MQ_bins(args)
    
    for i in range(len(_MAG_benchmark_blastn)):
        _MAG_benchmark_blastn[i] = parse_CRISPR_results(args,_MAG_benchmark_blastn[i],CONTIG2MAG)

    blastn_file = os.path.join('08_crisprcas','blastn','MAGS.all.virus.m6.gz')

    if os.path.exists(blastn_file):
        for i in range(len(_MAG_benchmark_blastn)):
            with gzip.open(blastn_file,'rt') as infile:
                for line in infile:
                    line = line.strip().split('\t')
                    phagebin = line[0]
                    phagebin = '_'.join(phagebin.split('_')[:2])

                    MAG_contig = line[1]
                    alignment = int(line[3])
                    MAG_contig_length = int(line[-1])
                    phagebin_length = int(line[-2])
                    
                    size_ratio = MAG_contig_length/phagebin_length
                    if alignment >= 500 and size_ratio > phagesize_ratios[i]:
                        MAG = CONTIG2MAG[MAG_contig]
                        phagemotherbin = phagebin.split('_')[1]

                        if not phagebin in medium_quality_bins:
                            continue

                        if not phagemotherbin in _MAG_benchmark_blastn[i][MAG]:
                            _MAG_benchmark_blastn[i][MAG][phagemotherbin] = {'spacer':None,'prophage':True,'both':None}
                        else:
                            if phagemotherbin in _MAG_benchmark_blastn[i][MAG] and _MAG_benchmark_blastn[i][MAG][phagemotherbin]['spacer'] is True:
                                _MAG_benchmark_blastn[i][MAG][phagemotherbin] = {'spacer':True,'prophage':True,'both':True}
    else:
        print('The blastn file:', blastn_file, 'does not exists!')
        sys.exit(1)

    ### Calculate Host-prediction accuracy     
    contig_host_results = {i:None for i in _MAG_benchmark_blastn}
    for i in _MAG_benchmark_blastn:
        contig_host_results[i] =  evaluate_host_consistency(_MAG_benchmark_blastn[i]) 

    host_prediction_results = return_host_consistency_table(contig_host_results,prediction_method='blastn')
    combined_table = np.array(host_prediction_results)
    combined_table = np.concatenate(combined_table,axis=0)
    fileout = os.path.join('08_crisprcas','crispr_prophage_results','blastn_host_predictions.tsv')
    np.savetxt(fileout, combined_table , delimiter='\t', header='', newline='\n', fmt='%s')


    for i in range(len(_MAG_benchmark_blastn)):

        ratioid = 'ratio_' + str(phagesize_ratios[i]) 
        outfile = os.path.join('08_crisprcas','crispr_prophage_results','MAG_crispr_prophage_summary.' + ratioid + '.blastn.txt')
        bacterial_genome_decision_contig = {MAG:{'spacer':None,'prophage':None,'both':None} for MAG in _MAG_benchmark_blastn[i]}
        with open(outfile,'w') as out:
            for MAG in _MAG_benchmark_blastn[i]:
                bacterial_cluster = MAG.split('_')[1]
                for phagemotherbin in _MAG_benchmark_blastn[i][MAG]:
                    MAG_sample = MAG.split('_')[0]
                    spacer = _MAG_benchmark_blastn[i][MAG][phagemotherbin]['spacer']
                    prophage =  _MAG_benchmark_blastn[i][MAG][phagemotherbin]['prophage']
                    both = _MAG_benchmark_blastn[i][MAG][phagemotherbin]['both']

                    if not str(phagemotherbin) in population_type:
                        continue
                    ctype = population_type[phagemotherbin]

                    if both is True: 
                        viral_association = 'both'
                    else:
                        if spacer is True:
                            viral_association = 'spacer'
                        else:
                            viral_association = 'prophage'


                    lineout = '\t'.join( [bacterial_cluster,MAG,phagemotherbin,ctype,viral_association] )
                    out.write(lineout+'\n')

                    ### 
                    if ctype in ['HQ-ref']:
                        bacterial_genome_decision_contig[MAG][viral_association] = True


        outfile = os.path.join('08_crisprcas','crispr_prophage_results','MAG_crispr_prophage_decision.' + ratioid + '.blastn.txt')
        with open(outfile,'w') as out:
            for MAG in bacterial_genome_decision_contig:
                bacterial_cluster = MAG.split('_')[1]
                viral_association = 'None'
                spacer = bacterial_genome_decision_contig[MAG]['spacer']
                prophage =  bacterial_genome_decision_contig[MAG]['prophage']
                both = bacterial_genome_decision_contig[MAG]['both']

                if prophage is True and spacer is True:
                    viral_association = 'both'
                    lineout = '\t'.join( [bacterial_cluster,MAG,viral_association] )
                    out.write(lineout+'\n')
                else:
                    if prophage is True:
                        viral_association = 'prophage'
                    elif spacer is True:
                        viral_association = 'spacer'
                    lineout = '\t'.join( [bacterial_cluster,MAG,viral_association] )
                    out.write(lineout+'\n')


if __name__ == "__main__":
    args = parser.parse_args()

    if not os.path.exists('08_crisprcas/crispr_prophage_results'):
        os.makedirs('08_crisprcas/crispr_prophage_results')

    ### Parsing results with alignments based on Blastn results
    parse_CRISPR_PROPHAGE_hits_blastn(args)

    ### Parsing results with alignments based on Fastani results
    parse_CRISPR_PROPHAGE_hits(args)




