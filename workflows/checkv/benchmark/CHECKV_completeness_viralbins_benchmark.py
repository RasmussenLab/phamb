#!/bin/python

import pandas as pd 
import numpy as np
import argparse
import os
import sys
import csv
from Bio import SeqIO
import gzip
from collections import defaultdict

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-c',help='VAMB directory')
parser.add_argument('-a',help='Annotation summary files directory')
parser.add_argument('-v', help='checkv directory')


def load_MGXMVX_HQ():


    def parse_fastani_hits(MGX_MVX_fastani_file,HQ_vOTUs):
        '''
        '''
        
        matching_MGX = dict()
        with open(MGX_MVX_fastani_file,'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                vOTU = os.path.basename(line[0]).replace('.fna','')
                MGXbin_name = os.path.basename(line[1]).replace('.fna','')

                if not vOTU in HQ_vOTUs:
                    continue
                ANI = float(line[2])
                percentage_of_vOTU = int(line[3])/int(line[4])*100

                if ANI >= 95 and percentage_of_vOTU >= 75:
                    MGX_bin = bin_class()
                    MGX_bin.vOTU = (vOTU,ANI,percentage_of_vOTU)
                    MGX_bin.contigs = dict()
                    matching_MGX[MGXbin_name] = MGX_bin

        return matching_MGX 

    def parse_cluster_file(clusters_file,MGXMVX_bins_clean):

        contig_to_bin = dict()
        with open(clusters_file,'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                contigname = line[1]
                clstr = line[0]
                binid = contigname.split('_')[0] + '_' + clstr
                if not binid in MGXMVX_bins_clean:
                    continue
                contig_to_bin[contigname] = binid

                MGXMVX_bins_clean[binid].contigs[contigname] = None

        return MGXMVX_bins_clean,contig_to_bin 


    def parse_blastn_hits(MGX_MVX_blastn_file,MGXMVX_bins_clean,contigs_overview):
        ''' 
        Record if all contigs in a cluster matches the vOTU assigned
        '''
        with open(MGX_MVX_blastn_file,'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                contigname, vOTU = line[0],line[1]
                if not contigname in contigs_overview:
                    continue

                if contigs_overview[contigname].vOTU == vOTU:
                    mgxbin = contigs_overview[contigname].mgxbin
                    MGXMVX_bins_clean[mgxbin].contigs[contigname] = vOTU
        return MGXMVX_bins_clean

    class contig_class:
        def __init__(self):
            pass
    
    class bin_class:
        def __init__(self):
            pass
        
    MGX_MVX_blastn_file = '06_blast_MGXMVX/MGX.2000.copsac.assemblies.vOTUs.new.m6'
    clusters_file = '05_binning/vamb_on_jgi/COPSAC/clusters.tsv'
    MGX_MVX_fastani_file = '06_fastani_MGXMVX/mvx_mgs.all.fastani.2.txt'

    ### Determine the HQ MVX genomes From COPSAC / Shiraz genomes
    p = 'combined_assemblies/checkv.out.new.vOTU/quality_summary.tsv'
    HQ_vOTUs = dict()
    for r in csv.DictReader(open(p), delimiter="\t"):
        vOTU = r["contig_id"]
        quality = r["checkv_quality"]
        completeness = r["completeness"]
        if quality in ['High-quality','Complete']:
            HQ_vOTUs[vOTU] = 'NA'

    MGXMVX_bins = parse_fastani_hits(MGX_MVX_fastani_file,HQ_vOTUs)
    
    ### Remame weird ending bin-names
    MGXMVX_bins_clean = dict()
    for binid in MGXMVX_bins:
        s = binid.split('_')
        if len(s) == 3:
            new_binid = '_'.join(s[:2])
        else:
            new_binid = binid
        MGXMVX_bins_clean[new_binid] =  MGXMVX_bins[binid]

    MGXMVX_bins_clean,contig_to_bin = parse_cluster_file(clusters_file,MGXMVX_bins_clean)

    contigs_overview = dict()
    for c in contig_to_bin:
        contig = contig_class()
        contig.mgxbin = contig_to_bin[c]
        contig.vOTU = MGXMVX_bins_clean[contig.mgxbin].vOTU[0]
        contigs_overview[c] = contig

    MGXMVX_bins_clean = parse_blastn_hits(MGX_MVX_blastn_file,MGXMVX_bins_clean,contigs_overview)

    ### Calculate how mnany contigs mapping to the vOTU in each MGX-bin
    for binid in MGXMVX_bins_clean:
        c = MGXMVX_bins_clean[binid].contigs
        nc = len(c)
        nc_mapping = sum([1 for i in c if not i is None])
        MGXMVX_bins_clean[binid].number_of_contigs = nc
        MGXMVX_bins_clean[binid].number_of_mapping_contigs = nc_mapping 

    return MGXMVX_bins_clean






def load_hmm_annotations(args):
    hmms = {}
    p = os.path.join(args.db, "hmm_db/genome_lengths.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        r["genomes"] = int(r["genomes"])
        r["cv"] = float(r["cv"])
        r["lengths"] = [int(_) for _ in r["lengths"].split(",")]
        hmms[r["hmm"]] = r
    return hmms


    
def read_checkv_refs(args):

    refs = {}
    p = os.path.join(args.db, "genome_db/checkv_reps.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        genome = Genome()
        genome.id = r["checkv_id"]
        genome.length = int(r["length"])
        genome.type = r["type"]
        genome.weight = (
            1 if r["type"] == "circular" else 2 if r["type"] == "genbank" else None
        )
        refs[genome.id] = genome
    return refs


def read_checkv_aai(args):
    
    aai_file = os.path.join(args.v,'tmp','aai.tsv')

    contig_aai_tophits = dict()
    with open(aai_file,'r') as infile:
        header = infile.readline().strip().split('\t')
        target_index = header.index('target')
        percent_index = header.index('percent_genes')
        percent_length_index = header.index('percent_length')
        aligned_genes_index = header.index('aligned_genes')
        aligned_length_index = header.index('aligned_length')
        identity_index = header.index('identity')
        score_index = header.index('score')
        for line in infile:
            line = line.strip().split('\t')
            contigname = line[0]
            target, aligned_genes, identity, score, aligned_length, aai_af = line[target_index],line[aligned_genes_index], line[identity_index], line[score_index], line[aligned_length_index], line[percent_length_index]

            if not contigname in contig_aai_tophits:
                contig_aai_tophits[contigname] = dict()
                contig_aai_tophits[contigname][target] = {'identity':identity,'aligned_length':aligned_length,'aai_af': aai_af,'score':score}
                largest_score = score
            else:
                ### Only contigs matching highly similar Genomes 
                if float(score)/float(largest_score) >= 0.75:
                    contig_aai_tophits[contigname][target] = {'identity':identity,'aligned_length':aligned_length,'aai_af': aai_af,'score':score}
    return contig_aai_tophits



def read_checkv_contamination(args):
    """[summary]

    Args:
        args ([type]): [description]

    Returns:
        [type]: [description]
    """

    contamination_file = os.path.join(args.v,'contamination.tsv')
    checkv_contigs = set()
    checkv_contig_contamination = dict()
    with open(contamination_file,'r') as infile:
        header = infile.readline().strip().split('\t')
        viral_genes_index = header.index('viral_genes')
        host_genes_index = header.index('host_genes')
        gene_count_index = header.index('total_genes')
        contiglen_index = header.index('contig_length')
        region_index = header.index('region_types')

        for line in infile:
            line = line.strip().split('\t')
            contigname = line[0]
            checkv_contigs.add(contigname)

            checkv_contig_contamination[contigname] = {'contiglen':line[contiglen_index],
            'total_genes':line[gene_count_index],
            'viral_genes':line[viral_genes_index],
            'host_genes':line[host_genes_index],
            'region_types':line[region_index]}

    return checkv_contig_contamination, checkv_contigs


def read_checkv_repeat(args):
    """[summary]

    Args:
        args ([type]): [description]

    Returns:
        [type]: [description]
    """

    contamination_file = os.path.join(args.v,'repeats.tsv')
    checkv_contig_repeat = dict()
    with open(contamination_file,'r') as infile:
        header = infile.readline().strip().split('\t')
        repeat_index = header.index('repeat_type')
        repeat_length_index = header.index('repeat_length')

        for line in infile:
            line = line.strip().split('\t')
            contigname = line[0]

            checkv_contig_repeat[contigname] = {'repeattype':line[repeat_index],
            'repeatlen':line[repeat_length_index]}
    return checkv_contig_repeat


def read_checkv_completeness(args):
    """[summary]

    Args:
        args ([type]): [description]

    Returns:
        [type]: [description]
    """

    contamination_file = os.path.join(args.v,'completeness.tsv')
    checkv_contig_completeness = dict()
    with open(contamination_file,'r') as infile:
        header = infile.readline().strip().split('\t')
        tophit_aai_index = header.index('aai_id')
        tophit_af_index = header.index('aai_af')
        tophit_index = header.index('aai_top_hit')
        aai_conf_index = header.index('aai_confidence')
        aai_error_index = header.index('aai_error')
        viral_length_index = header.index('viral_length')
        completeness_index = header.index('aai_completeness')
        reference_length_index = header.index('aai_expected_length')

        for line in infile:
            line = line.strip().split('\t')
            contigname = line[0]

            checkv_contig_completeness[contigname] = {'virallen':line[viral_length_index],
            'tophit_aai':line[tophit_aai_index],
            'tophit_af':line[tophit_af_index],
            'tophit':line[tophit_index],
            'tophit_size':line[reference_length_index],
            'aai_conf':line[aai_conf_index],
            'aai_error':line[aai_error_index]
            }

    return checkv_contig_completeness



def parse_VAMB_clusters(args,checkv_contigs):

    clusterfile = os.path.join(args.c,'clusters.tsv')
    bins = dict()
    with open(clusterfile,'r') as infile:
        for line in infile:
            cluster, contig = line.strip().split()

            if not contig in checkv_contigs:
                continue
            sample = contig.split('_')[0]
            binid = sample + '_' + cluster
            
            if binid not in bins:     
                mothercluster = binid.split('_')[1]
                binobject = VAMB_bin_class(binid,mothercluster)
                binobject.contigs.append(contig)
                bins[binid] = binobject
            else:
                binobject = bins[binid]
                binobject.contigs.append(contig)

    return bins



def flush_class(VAMB_bin):
    VAMB_bin.checkv_ref = []
    VAMB_bin.binsize = 0
    VAMB_bin.contig_lengths = []
    VAMB_bin.checkv_ref = []
    VAMB_bin.checkv_ref_size = dict()
    VAMB_bin.host_type = []
    VAMB_bin.confidence = []
    VAMB_bin.aai = []
    VAMB_bin.aai_error = []
    VAMB_bin.aai_af= []
    VAMB_bin.circ = []
    VAMB_bin.circ_len = []

    return VAMB_bin


class VAMB_bin_class:
     def __init__(self,bin_name,mothercluster):
        self.bin = bin_name
        self.binsize = 0
        self.mothercluster = mothercluster
        self.contigs = []
        self.contig_lengths = []
        self.checkv_ref = []
        self.checkv_ref_size = dict()
        self.host_type = []
        self.confidence = []
        self.aai = []
        self.aai_error = []
        self.aai_af= []
        self.circ = []
        self.circ_len = []
        self.hmms = []
    

def get_confidence(error):
    confidence = (
            "low"
            if error == "NA"
            else "high"
            if error <= 5
            else "medium"
            if error <= 10
            else "low"
        )
    return confidence

def get_cluster_splits(VAMB_bin, contig_aai_tophits,avg_aai_cutoff=80, ref_coverage_cutoff = 50):
    references = set(VAMB_bin.checkv_ref)
    references = [r for r in references if r != 'NA']
    cluster_splits = dict()
    j = 0
    for checkv_ref in references:
        nearest_ref_size = round(float(VAMB_bin.checkv_ref_size[checkv_ref]),2)
        genome_size_in_bp = 0
        aligned_length = 0
        aai = []
        aai_af = []
        aai_error = []
        indices = []
        for i,c in enumerate(VAMB_bin.contigs):

            if not c in contig_aai_tophits: # Probably noise in these cases
                continue

            if VAMB_bin.checkv_ref[i] == checkv_ref:
                aligned_length += int(contig_aai_tophits[c][checkv_ref]['aligned_length'])
                genome_size_in_bp += int(VAMB_bin.contig_lengths[i])
                aai += [float(contig_aai_tophits[c][checkv_ref]['identity'])]
                aai_af += [float(contig_aai_tophits[c][checkv_ref]['aai_af'])]

                if VAMB_bin.aai_error[i] == 'NA':
                    aai_error += [np.nan]
                else:
                    aai_error += [float(VAMB_bin.aai_error[i])]
                indices += [i]

        avg_aai = round(np.mean(aai),2)
        avg_aai_error = round(np.nanmean(aai_error),2)
        confidence = get_confidence(avg_aai_error)
        aai_af = round(np.mean(aai_af),2)
        aai_completeness = round( (genome_size_in_bp/nearest_ref_size)*100,2)
        
        if avg_aai >= avg_aai_cutoff and aai_af >= ref_coverage_cutoff:
            j += 1
            cluster_splits[VAMB_bin.bin + '_' + str(j)] = {'confidence':confidence,
            'avg_aai_error':avg_aai_error,
            'contig_indices':indices,
            'avg_aai':avg_aai,
            'aai_af':aai_af,
            'genome_coverage_in_bp':genome_size_in_bp,
            'aai_completeness':aai_completeness,
            'nearest_ref_size':nearest_ref_size,
            'checkv_ref':checkv_ref}
    return cluster_splits




def viral_hmm_annotation(VAMB_bins,annotation_path,hmms):

    contig_to_bin  = dict()
    for VAMB_bin in VAMB_bins:
        for c in VAMB_bins[VAMB_bin].contigs:
            contig_to_bin[c] = VAMB_bin

    for r in csv.DictReader(open(annotation_path), delimiter="\t"):
        if r["target_hmm"] in hmms and hmms[r["target_hmm"]]["genomes"] >= 10:
            c = r["contig_id"]
            vb = contig_to_bin[c]
            VAMB_bins[vb].hmms.append(r["target_hmm"])
    return VAMB_bins

def viral_hmm_calculation(VAMB_bin,hmms):

    # estimate minimum completeness
        # check at least 1 hmm
    if len(VAMB_bin.hmms) == 0:
        return None
    else:
        # fetch genome length (viral region only)
        query_length = VAMB_bin.binsize

        # get weighted completeness range
        comps = []
        for hmm in VAMB_bin.hmms:
            len_q1 = np.quantile(hmms[hmm]["lengths"], 0.05)
            len_q2 = np.quantile(hmms[hmm]["lengths"], 0.95)
            comp_q1 = round(min([100 * query_length / len_q2, 100.0]), 2)
            comp_q2 = round(min([100 * query_length / len_q1, 100.0]), 2)
            cv = hmms[hmm]["cv"]
            weight = min([1 / cv, 50]) if cv != 0 else 50
            comps.append([weight, comp_q1, comp_q2])
        weight_total = sum([weight for weight, comp_q1, comp_q2 in comps])
        comp_lower = (
            sum([weight * comp_q1 for weight, comp_q1, comp_q2 in comps]) / weight_total
        )
        comp_upper = (
            sum([weight * comp_q2 for weight, comp_q1, comp_q2 in comps]) / weight_total
        )
        hmm_completeness = (round(comp_lower,2), round(comp_upper,2))
        return hmm_completeness


def identify_best_representative(VAMB_bin,contig_aai_tophits,refs):

    references = set(VAMB_bin.checkv_ref)
    references = [r for r in references if r != 'NA']

    cluster_splits = dict()
    j = 0
    best_score, best_ref = 0, None
    for checkv_ref in references:
        nearest_ref_size = round(float(VAMB_bin.checkv_ref_size[checkv_ref]),2)
        genome_size_in_bp = 0
        aligned_length = 0
        aai = []
        aai_af = []
        aai_error = []
        indices = []
        best_target = []
        for i, c in enumerate(VAMB_bin.contigs):

            if not c in contig_aai_tophits: # Probably noise in these cases
                continue

            if checkv_ref in contig_aai_tophits[c]:
                aligned_length += int(contig_aai_tophits[c][checkv_ref]['aligned_length'])
                genome_size_in_bp += int(VAMB_bin.contig_lengths[i])
                aai += [float(contig_aai_tophits[c][checkv_ref]['identity'])]
                aai_af += [float(contig_aai_tophits[c][checkv_ref]['aai_af'])]
                if VAMB_bin.aai_error[i] == 'NA':
                    aai_error += [np.nan]
                else:
                    aai_error += [float(VAMB_bin.aai_error[i])]
                indices += [i]
        
        avg_aai = round(np.mean(aai),2)
        avg_aai_error = round(np.nanmean(aai_error),2)
        confidence = get_confidence(avg_aai_error)
        aai_af = round(np.mean(aai_af),2)
        aai_completeness = round( (genome_size_in_bp/nearest_ref_size)*100,2)
        score = round(avg_aai * aligned_length / 100,2)
        if score > best_score:
            best_score = score
            best_ref = checkv_ref
        j += 1
        cluster_splits[checkv_ref] = {'confidence':confidence,
        'avg_aai_error':avg_aai_error,
        'contig_indices':indices,
        'avg_aai':avg_aai,
        'aai_af':aai_af,
        'genome_coverage_in_bp':genome_size_in_bp,
        'nearest_ref_size':nearest_ref_size,
        'aai_completeness':aai_completeness,
        'checkv_ref':checkv_ref,
        'score':score}

    ### Estimat New expected lengths based on all contigs in the highest scoring Contig-combination
    # If one exists.
    if not best_ref is None:
        indices = cluster_splits[best_ref]['contig_indices']
        genome_size_in_bp = cluster_splits[best_ref]['genome_coverage_in_bp']

        all_top_targets = []
        all_scores = []
        for i in indices:
            c = VAMB_bin.contigs[i]
            if not c in contig_aai_tophits: # Probably noise in these cases
                continue
            for r in contig_aai_tophits[c]:
                all_top_targets += [r]
                all_scores += [float(contig_aai_tophits[c][r]['score'])]
        lengths = [refs[r].length for r in all_top_targets]
        weights = [refs[r].weight * all_scores[i] for i,r in enumerate(all_top_targets)]
        new_expected_length = sum(l * w for l, w in zip(lengths, weights)) / sum(weights)

        aai_completeness = round( (genome_size_in_bp/new_expected_length)*100,2)
        cluster_splits[best_ref]['aai_completeness'] = aai_completeness
        cluster_splits[best_ref]['nearest_ref_size'] =  new_expected_length

        ### Note the remainnder/non-mapping contigs 
        unmapped_contig = [ float(VAMB_bin.contig_lengths[i]) for i,c in enumerate(VAMB_bin.contigs) if i not in indices ]
        bp_unmapped_contig = sum(unmapped_contig)
        cluster_splits[best_ref]['bp_unmapped_contig'] = bp_unmapped_contig
        

    return cluster_splits, best_ref


def gather_contig_annotation(VAMB_bin,checkv_results):
    VAMB_bin = flush_class(VAMB_bin)
    for c in VAMB_bin.contigs:
        conf = checkv_contig_completeness[c]['aai_conf']
        aai = checkv_contig_completeness[c]['tophit_aai']
        aai_error = checkv_contig_completeness[c]['aai_error']
        aai_af = checkv_contig_completeness[c]['tophit_af']
        ref = checkv_contig_completeness[c]['tophit']
        ref_size = checkv_contig_completeness[c]['tophit_size']
        host = checkv_contig_contamination[c]['region_types']
        contig_length = checkv_contig_contamination[c]['contiglen']
        viral_length = checkv_contig_completeness[c]['virallen']
        if viral_length != 'NA':
            contig_length = viral_length
        circ, circ_len = None, None
        if c in checkv_contig_repeat:
            circ = checkv_contig_repeat[c]['repeattype']
            circ_len = checkv_contig_repeat[c]['repeatlen']
        
        ### 
        VAMB_bin.binsize += float(contig_length)
        VAMB_bin.checkv_ref_size[ref] = ref_size
        VAMB_bin.checkv_ref += [ref]
        VAMB_bin.contig_lengths += [contig_length]
        VAMB_bin.confidence += [conf]
        VAMB_bin.aai += [aai]
        VAMB_bin.aai_error += [aai_error]
        VAMB_bin.aai_af += [aai_af]
        VAMB_bin.circ += [circ]
        VAMB_bin.circ_len += [circ_len]
        VAMB_bin.host_type += [host]
    
    return VAMB_bin


def return_circular_splits(VAMB_bin):

    circular_indices = []
    for i, c in enumerate(VAMB_bin.contigs):
        if VAMB_bin.circ[i] != 'NA':
            circular_indices.append(i)

    cluster_splits = dict()
    j = 1 
    for i in circular_indices:
        checkv_ref = VAMB_bin.checkv_ref[i]

        if VAMB_bin.checkv_ref_size[checkv_ref] == 'NA':
            aai_af, aai_completeness, aai_error, ref_size, aai,confidence = 'NA','NA','NA','NA','NA','NA'
        else:
            ref_size = round(float(VAMB_bin.checkv_ref_size[checkv_ref]),2)
            aai_error = round(float(VAMB_bin.aai_error[i]),2)
            aai_af = round( float(VAMB_bin.aai_af[i]),2)
            aai_completeness = round( (float( VAMB_bin.contig_lengths[i] )/ ref_size )*100 ,2)
            aai = round(float(VAMB_bin.aai[i]),2)
            confidence = VAMB_bin.confidence[i]

        cluster_splits[VAMB_bin.bin + '_' + str(j) + '_' + 'circular'] = {'confidence':confidence,
        'avg_aai_error':aai_error,
        'contig_indices':[i],
        'avg_aai': aai,
        'aai_af':  aai_af ,
        'genome_coverage_in_bp':VAMB_bin.contig_lengths[i],
        'nearest_ref_size':ref_size,
        'aai_completeness': aai_completeness,
        'checkv_ref':checkv_ref,
        'bp_unmapped_contig':0}
        j +=1

    return cluster_splits

def return_single_contig_split(VAMB_bin):
    checkv_ref = VAMB_bin.checkv_ref[0]

    if VAMB_bin.checkv_ref_size[checkv_ref] == 'NA':
        aai_af, aai_completeness, aai_error, ref_size, aai,confidence = 'NA','NA','NA','NA','NA','NA'
    else:
        
        ref_size = round(float(VAMB_bin.checkv_ref_size[checkv_ref]),2)
        aai_af = aai_af = round( float(VAMB_bin.aai_af[0]),2)
        aai_completeness = round( (float( VAMB_bin.contig_lengths[0] ) / ref_size )*100 ,2)
        aai = round(float(VAMB_bin.aai[0]),2)
        aai_error = VAMB_bin.aai_error[0]
        confidence = VAMB_bin.confidence[0]
        if not aai_error == 'NA':
            aai_error = round(float(VAMB_bin.aai_error[0]),2)

    
    cluster_splits = dict()
    cluster_splits[VAMB_bin.bin + '_' +str(1)] = {'confidence':confidence,
    'avg_aai_error':aai_error,
    'contig_indices':[0],
    'avg_aai': aai,
    'aai_af':  aai_af ,
    'genome_coverage_in_bp':VAMB_bin.contig_lengths[0],
    'nearest_ref_size':ref_size,
    'aai_completeness': aai_completeness,
    'checkv_ref':checkv_ref,
    'bp_unmapped_contig':0}

    return cluster_splits



def evaluate_VAMB_cluster(VAMB_bin,checkv_results,refs,hmms):
    """[summary]

    Args:
        VAMB_bin ([type]): [description]
        checkv_results ([type]): [description]

    Returns:
        [type]: [description]
    """

    contig_aai_tophits = checkv_results['contig_aai_tophits']

    ### Parse contig information
    VAMB_bin = gather_contig_annotation(VAMB_bin,checkv_results)

    ### Check HMM completeness

    hmm_completeness = viral_hmm_calculation(VAMB_bin,hmms)

    ### Check for Circular Viral genomes
    circular_splits = None
    if 'DTR' in VAMB_bin.circ or 'ITR' in VAMB_bin.circ:
        circular_splits = return_circular_splits(VAMB_bin)

    ### Single contig Bin - result will be what CheckV decided.
    if len(VAMB_bin.contigs) == 1:
        cluster_splits = return_single_contig_split(VAMB_bin)

    else:

        ### First check for bins containing mix of Genomes
        #cluster_splits = get_cluster_splits(VAMB_bin, contig_aai_tophits,avg_aai_cutoff=90, ref_coverage_cutoff=50)
        cluster_splits = dict()
        ### If no splits exists, it might be a case with contigs matching highly similar references but actually belong together.
        # In this, case we can do a majority count to find the best complete representative genome possible
        if len(cluster_splits) == 0:
            best_cluster_splits,best_ref = identify_best_representative(VAMB_bin,contig_aai_tophits,refs)
            
            if not best_ref is None:
                cluster_splits = dict()
                cluster_splits[VAMB_bin.bin + '_' +str(1)] = best_cluster_splits[best_ref]
            else:
                cluster_splits = {}
    
    ### Merge circular splits with the rest
    if not circular_splits is None:
        cluster_splits = {**circular_splits, **cluster_splits}
    return cluster_splits, hmm_completeness


def write_completeness_file(cluster_splits,VAMB_bin,hmm_completeness,out, HQ_MGXMVX_bins_clean):
    '''
    1. Write out completeness information similar to CheckV 
    2. [Special for Benchmark] - Note how many of the Contigs matching the vOTU are part of the refined VAMB-bin
    '''
    if hmm_completeness is None:
        hmm_completeness = ('NA','NA')

    if not cluster_splits is None: 

        for split in cluster_splits:
            
            contigs = VAMB_bin.contigs

            contigs_in_vOTU = HQ_MGXMVX_bins_clean[VAMB_bin.bin].contigs
            contigs_in_vOTU = [c for c in contigs_in_vOTU if not contigs_in_vOTU[c] is None]
            contigs_included = [contigs[i] for i in cluster_splits[split]['contig_indices'] ]
            vOTU_contigs_in_refined_bin = sum([1 for c in contigs_included if c in contigs_in_vOTU])
            
            ### Total number of contigs, contigs matching the vOTU and contigs in the refined bin matching the vOTU
            
            contigs_included = ';'.join(contigs_included)
            lineout = [VAMB_bin.bin,
            split, 
            cluster_splits[split]['confidence'], 
            cluster_splits[split]['avg_aai_error'], 
            len(cluster_splits[split]['contig_indices']),
            cluster_splits[split]['avg_aai'],
            cluster_splits[split]['aai_af'],
            cluster_splits[split]['genome_coverage_in_bp'],
            cluster_splits[split]['nearest_ref_size'],
            cluster_splits[split]['aai_completeness'],
            cluster_splits[split]['checkv_ref'],
            cluster_splits[split]['bp_unmapped_contig'],
            hmm_completeness[0], #Lower bound HMM estimate
            hmm_completeness[1], #Upper bound HMM estimate
            vOTU_contigs_in_refined_bin,
            len(contigs_in_vOTU),
            len(contigs),
            contigs_included]
            lineout = [str(x) for x in lineout]
            out.write('\t'.join(lineout)+'\n')
    #else:
       # lineout = [binid+'_'+'NA'] + ['NA' for i in range(14)]
       # out.write('\t'.join(lineout)+'\n')


class Args:
    def __init__(self,v,c,o,db):
        self.v = v 
        self.o = o
        self.c = c
        self.db = db


args = Args('07_binannotation/checkv/VAMB_contigs','05_binning/vamb_on_jgi/COPSAC','lel','/home/projects/cpr_10006/projects/phamb/databases/checkv/checkv-db-v0.6')
class Genome:
    def __init__(self):
        pass



if __name__ == "__main__":
    args = parser.parse_args()

    refs = read_checkv_refs(args)
    hmms = load_hmm_annotations(args) # HMM database
    annotation_path = os.path.join(args.v,"tmp", "gene_annotations.tsv")

    ### Define the overlapping MGX bins 
    # Then determine how many of the individual contigs actually match their according MVX genome
    HQ_MGXMVX_bins_clean = load_MGXMVX_HQ()

    checkv_contig_contamination, checkv_contigs = read_checkv_contamination(args)
    checkv_contig_repeat = read_checkv_repeat(args)
    checkv_contig_completeness = read_checkv_completeness(args)
    contig_aai_tophits = read_checkv_aai(args)

    VAMB_bins = parse_VAMB_clusters(args, checkv_contigs)

    ### Attach HMM hits to each VAMB bin
    VAMB_bins = viral_hmm_annotation(VAMB_bins,annotation_path,hmms)

    checkv_results = {'repeats':checkv_contig_repeat,'completeness':checkv_contig_completeness,'contamination':checkv_contig_contamination,'contig_aai_tophits':contig_aai_tophits}

    ### Write results
    fileout = os.path.join(args.v,'VAMB_completeness_onlyHQMVX.AF75.tsv')
    with open(fileout,'w') as out:

        header = ['vamb_bin',
        'vamb_split',
        'aai_confidence',
        'aai_error',
        'ncontigs',
        'aai_id',
        'aai_af',
        'viral_length',
        'aai_expected_length',
        'aai_completeness',
        'aai_top_hit',
        'unmapped_bp',
        'hmm_lower_completeness',
        'hmm_upper_completeness',
        'vOTU_contigs_after_refine',
        'vOTU_contigs_before_refine',
        'contigs_total_in_MGX_bin',
        'contigs_in_genome']
        out.write('\t'.join(header)+'\n')

        for binid in list(HQ_MGXMVX_bins_clean):

            VAMB_bin = VAMB_bins[binid]
            contigs = VAMB_bin.contigs

            cluster_splits,hmm_completeness = evaluate_VAMB_cluster( VAMB_bin, checkv_results ,refs,hmms)
            write_completeness_file(cluster_splits,VAMB_bin,hmm_completeness,out,HQ_MGXMVX_bins_clean)