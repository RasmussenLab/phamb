#!/bin/Python
### Python code to parse .npz into depth matrix in .txt format with contignames and contiglengths
import numpy as np
import sys 
import argparse
import os
import scipy.stats
import multiprocessing as mp
import gzip
import pickle

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
   ''')
parser.add_argument('-d', help='JGI depth matrix file')
parser.add_argument('-l', help='List of JGI coverage files')
parser.add_argument('-m', help='checkm-mags file')
parser.add_argument('-c', help='Cluster file - will be used to extract specific contigs of interest',default=None)
parser.add_argument('-o', help='file outbase')



class bin_abundance:
     def __init__(self,bin_name,mothercluster):
        self.bin = bin_name
        self.binsize = 0
        self.mothercluster = mothercluster
        self.contigs = set()
        self.binRPKM = None
        self.binRPM = None



def parse_prophages(blastn_file,population_type):
    
    prophages = dict()
    with gzip.open(blastn_file,'rt') as infile:
        for line in infile:
            line = line.strip().split('\t')
            phagebin, MAG_contig, sequence_identity, alignment = line[:4]
            phagesize, MAG_contig_size = line[-2],line[-1]
            phagecoverage = int(alignment)/int(phagesize) 

            
            
            if phagecoverage >= 0.80 and int(alignment) >= 4000:
                phagemotherbin = phagebin.split('_')[1]

                if not phagemotherbin in population_type:
                    continue

                ctype = population_type[phagemotherbin]
                if ctype in ['HQ-ref','Grey-matter']:
                    if len(phagebin.split('_')) == 3:
                        phagebin = '_'.join(phagebin.split('_')[:2])
                    if not MAG_contig in prophages:
                        prophages[MAG_contig] = dict()
                        prophages[MAG_contig][phagemotherbin] = (phagebin,int(alignment),float(sequence_identity),phagecoverage,ctype)
                    else:
                        if phagemotherbin in prophages[MAG_contig]:
                            current_alignment = prophages[MAG_contig][phagemotherbin][1]
                            if int(alignment) > current_alignment:
                                prophages[MAG_contig][phagemotherbin] = (phagebin,int(alignment),float(sequence_identity),phagecoverage,ctype)
                        else:
                            prophages[MAG_contig][phagemotherbin] = (phagebin,int(alignment),float(sequence_identity),phagecoverage,ctype)
    return prophages

def calculate_avg_motherbinRPKM(bin_abundances,motherbins):
        
    ### Calculate a combined mothercluster abundance
    motherbins_abundances = {}
    for clustername in motherbins:
        matrixsubset = []
        bins = motherbins[clustername]
        for binid in bins:
            matrixsubset.append(bin_abundances[binid].binRPKM)
        matrixsubset = np.array(matrixsubset)
        
        ### Avg. binabundances 
        clusterabundances = np.nanmean(matrixsubset, axis=0)
        clusterabundances[np.isnan(clusterabundances)] = 0
        motherbins_abundances[clustername] = clusterabundances
    return motherbins_abundances

def calculate_avg_motherbinRPM(bin_abundances,motherbins):
        
    ### Calculate a combined mothercluster abundance
    motherbins_abundances = {}
    for clustername in motherbins:
        matrixsubset = []
        bins = motherbins[clustername]
        for binid in bins:
            matrixsubset.append(bin_abundances[binid].binRPM)
        matrixsubset = np.array(matrixsubset)
        
        ### Avg. binabundances 
        clusterabundances = np.nanmean(matrixsubset, axis=0)
        clusterabundances[np.isnan(clusterabundances)] = 0
        motherbins_abundances[clustername] = clusterabundances
    return motherbins_abundances

def get_bin_abundances(args,sample_list, breadth_of_coverage, contig_order):

    ### parse contigs of cluster file 
    clusterfile = args.c
    magfile = args.m

    breadth_of_coverage = pickle.load( open('10_abundances/contig_breadth_of_coverage.pickle', "rb" ) )

    clusterfile = '05_binning/vamb_on_jgi_v3/HMP2/clusters.tsv'
    magfile = '07_binannotation/bacteria/checkm/HMP2.all.MQNC.MAGS'

    contigfile = '05_binning/vamb_on_jgi_v3/HMP2/contigs.npz'
    contig_order = np.load(contigfile)['arr_0']

    contiglengths =  '05_binning/vamb_on_jgi_v3/HMP2/lengths.npz'
    contig_lengths = np.load(contiglengths)['arr_0']

    contig_lengths_lookup = {}
    matrix_indices_dict = {}
    for i in range(len(contig_order)):
        contigname = contig_order[i]
        contiglength = contig_lengths[i]
        matrix_indices_dict[contigname] = i
        contig_lengths_lookup[contigname] = contiglength
    
    ### Load MAGS
    mags = set()
    with open(magfile,'r') as infile:
        infile.readline()
        for line in infile:
            MAG = line.strip().split('\t')[0]
            mags.add(MAG.split('_')[1])

    contig_to_mag = dict()
    with open(clusterfile,'r') as infile:
        for line in infile:
            line = line.strip().split()
            mothercluster, contig = line[0],line[1]
            binid = contig.split('_')[0] + '_' + mothercluster
            if not mothercluster in mags:
                continue
            contig_to_mag[contig] = binid


    population_type = dict()
    with open('07_binannotation/checkv/VAMB_bins/population_types.tsv','r') as infile:
        for line in infile:
            phagemotherbin, ctype = line.strip().split()
            population_type[phagemotherbin] = ctype

    ### Load Integrated Prophage contigs 
    blastn_file = '13_novoclusters/blastn/MAGS.all.virus.m6.gz'
    prophage_contigs = parse_prophages(blastn_file,population_type)

    ### Determine which viral motherbins we need to calculate avg. RPKM for. 
    prophage_motherbins = []
    for c in prophage_contigs:
        mbins = [prophage_motherbins.append(x) for x in prophage_contigs[c] ]
    prophage_motherbins = set(prophage_motherbins)
    prophage_MAGs = [contig_to_mag[c] for c in prophage_contigs.keys()]

    bin_abundances = dict()
    with open(clusterfile,'r') as infile:
        for line in infile:
            line = line.strip().split()
            mothercluster, contig = line[0],line[1]
            binid = contig.split('_')[0] + '_' + mothercluster
            if not binid in prophage_MAGs and not mothercluster in prophage_motherbins:
                continue
            contiglen = contig_lengths_lookup[contig]
            if not binid in bin_abundances:     
                mothercluster = binid.split('_')[1]
                binobject = bin_abundance(binid,mothercluster)
                binobject.contigs.add(contig)
                binobject.binsize += contiglen
                bin_abundances[binid] = binobject
            else:
                binobject = bin_abundances[binid]
                binobject.contigs.add(contig)
                binobject.binsize += contiglen
    
    ### Parse depth matrix and calculate weighted mean
    matrixfile = args.d
    matrixfile = '05_binning/vamb_on_jgi_v3/HMP2/RPKM.npz'
    matrix = np.load(matrixfile)['arr_0']


    samplefile = '05_binning/vamb_on_jgi_v3/HMP2/samples.npz'
    sample_order = np.load(samplefile)['arr_0']
    sample_order_lookup = {k:i for i,k in enumerate(sample_order) }
    sample_order_lookup_indices = {i:k for i,k in enumerate(sample_order) }    
    sample_indices = np.array( [sample_order_lookup[x[1]] for x in sample_list]) 
    new_sample_order = [sample_order_lookup_indices[i] for i in sample_indices]
    nsamples = len(sample_indices)
    matrix_subset = matrix[:,sample_indices]

    ### Calculate Bin-abundance across all samples - if coverage of a contig is too low, set depth to zero.
    for binid in bin_abundances:

        binobject = bin_abundances[binid]
        bin_contigs = binobject.contigs

        contig_indices = [ matrix_indices_dict[contigname] for contigname in bin_contigs]
        tmp_binabundances = matrix_subset[np.array(contig_indices),]
        
        ### Make sure to set the RPKM to zero in samples with low coverage for each contig, 
        # If the contig breadth_of_coverage is > 75% set to False and keep the RPKM value for the contig. Else set to Zero.
        contig_indices = []
        for j,contig in enumerate(bin_contigs):
            ones_boolean = np.ones((1, nsamples) , dtype = bool)[0]
            for i,ix in enumerate(sample_indices):
                s = sample_order_lookup_indices[ix]
                if contig in breadth_of_coverage[s]:
                    ones_boolean[i] = False
            tmp_binabundances[j,ones_boolean] = 0

        #bin_weighted_sample_averages = np.average(binabundances,axis=0, weights=contig_lengths)
        #bin_sample_sums = np.sum(binabundances,axis=0)
        
        tmp_binabundances[tmp_binabundances==0] = np.nan
        bin_sample_RPKM = np.nanmean(tmp_binabundances,axis=0)
        bin_sample_RPM = bin_sample_RPKM*(binobject.binsize/1000)
        binobject.binRPKM = bin_sample_RPKM
        binobject.binRPM = bin_sample_RPM
    
    ### Organise bins belonging to each Mothercluster
    motherbins = dict()
    for binid in bin_abundances:
        clustername = bin_abundances[binid].mothercluster
        if not clustername in motherbins:
            motherbins[clustername] = set([binid])
        else:
            motherbins[clustername].add(binid)
    
    motherbins_avgRPKM = calculate_avg_motherbinRPKM(bin_abundances,motherbins)
    motherbins_avgRPM = calculate_avg_motherbinRPM(bin_abundances,motherbins)

    ### Remember bin_abundances is a dict-object of class-objects
    ### motherbins_abundances is a standard dict 

    ### Now parse Through Prophage contigs - calculate the RPKM ratio of Free-Phage vs. contig w. phage integrated in that sample.
    prophage_RPKM_table = []
    prophage_samples = dict()
    for contig in prophage_contigs:
        s = contig.split('_')[0]

        if not s in new_sample_order:
            continue
        six = new_sample_order.index(s)
        cix = matrix_indices_dict[contig]

        contig_RPKM = matrix_subset[cix, six]
        MAG = contig_to_mag[contig]
        MAG_motherbin = MAG.split('_')[1]
        MAG_RPKM = motherbins_avgRPKM[MAG_motherbin][six]

        prophage_ratios[contig] = dict()
        for phagemotherbin in prophage_contigs[contig]:

            ### 
            phagebin = prophage_contigs[contig][phagemotherbin][0]
            #phage_RPKM = motherbins_avgRPKM[phagemotherbin][six]
            phage_RPKM = motherbins_avgRPKM[phagemotherbin][six]
            #bin_abundances[phagebin].binRPKM[six]

            if phage_RPKM >0:
                cov = prophage_contigs[contig][phagemotherbin][3]
                contig_phage_ratio = contig_RPKM/phage_RPKM
                MAG_phage_ratio = MAG_RPKM/phage_RPKM
                prophage_ratios[contig][phagebin] = (contig_phage_ratio,cov,MAG)
                row = [contig,s, MAG, phagebin,contig_RPKM,MAG_RPKM,phage_RPKM,phagemotherbin]
                prophage_RPKM_table += [row]
   
                if not phagemotherbin in prophage_samples:
                    prophage_samples[phagemotherbin] = [(s,MAG)]
                else:
                    prophage_samples[phagemotherbin] += [(s,MAG)]                  

    fileout = os.path.join('10_abundances','prophage_contig_RPKM.txt')

    with open(fileout,'w') as out:
        for row in prophage_RPKM_table:
            row = [str(x) for x in row]
            line = '\t'.join(row)
            out.write(line+'\n')

    ### Now we need to calculate what the corresponding 'Free' Phage
    # We basically know that the  
    free_pratio = []

    free_phage_RPKM_table = []
    #freephage_ratios = dict()
    for phagemotherbin in prophage_samples:
        hosts = set([p[1].split('_')[1] for p in prophage_samples[phagemotherbin]])
        pro_phage_samples = [p[0] for p in prophage_samples[phagemotherbin]]
        non_prophage_samples = [s for s in new_sample_order if not s in pro_phage_samples]

        for sample in non_prophage_samples:
            six = new_sample_order.index(sample)
            phagebin = sample +'_'+phagemotherbin 

            phage_RPKM = motherbins_avgRPKM[phagemotherbin][six]

            if not phage_RPKM >0:
                continue

            for MAG_motherbin in hosts:
                MAG = sample + '_' + MAG_motherbin

                MAG_RPKM = motherbins_avgRPKM[MAG_motherbin][six]

                if not MAG_RPKM >0 :
                    continue

                row = [sample, MAG, phagebin,MAG_RPKM,phage_RPKM,phagemotherbin]
                free_phage_RPKM_table += [row]

    fileout = os.path.join('10_abundances','freephage_contig_RPKM.txt')
    with open(fileout,'w') as out:
        for row in free_phage_RPKM_table:
            row = [str(x) for x in row]
            line = '\t'.join(row)
            out.write(line+'\n')



    return bin_abundances, motherbins_avgRPKM, motherbins_avgRPM



def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h



def removekey(d, key):
        r = dict(d)
        del r[key]
        return r


def return_RPKM_matrix(motherbins_avgRPKM,motherbinorder):

    rpkmmatrix = []
    clusters = []
    for tup in motherbinorder:
        cluster, val = tup
        clusters.append(cluster)
        rpkmmatrix.append(motherbins_avgRPKM[cluster])
    rpkmmatrix = np.array(rpkmmatrix)
    rpkmmatrix = rpkmmatrix.astype('float32').round(2)

    clusters = np.array(clusters).reshape(-1,1)
    rpkmmatrix = np.column_stack((clusters,rpkmmatrix))
    return rpkmmatrix



def write_matrices(args ,samples , motherbins_avgRPKM, motherbins_avgRPM):
    
    outbase = args.o

    mbin_median = [ (cluster,np.mean(motherbins_avgRPKM[cluster])) for cluster in motherbins_avgRPKM]
    mbin_median = sorted(mbin_median , key= lambda x: x[1],reverse=True)

    ### Prevalence 
    mbin_prevalence =  [ (cluster,np.sum(motherbins_avgRPKM[cluster] >0) ) for cluster in motherbins_avgRPKM]
    mbin_prevalence = sorted(mbin_prevalence , key= lambda x: x[1],reverse=True)

    rpkmmatrix = return_RPKM_matrix(motherbins_avgRPKM, mbin_median)
    header = ['mothercluster'] + list(samples[0][1:])
    header = '\t'.join(header)
    fileout1 = os.path.join(outbase,'VAMBV3.MAG.breadthcov.abundances.RPKM.txt')
    np.savetxt(fileout1, rpkmmatrix , delimiter='\t', header=header, newline='\n', fmt='%s')

    rpmmatrix = return_RPKM_matrix(motherbins_avgRPM, mbin_median)
    header = ['mothercluster'] + list(samples[0][1:])
    header = '\t'.join(header)
    fileout1 = os.path.join(outbase,'VAMBV3.MAG.breadthcov.abundances.RPM.txt')
    np.savetxt(fileout1, rpmmatrix , delimiter='\t', header=header, newline='\n', fmt='%s')






def read_write_MAG_taxonomy(args):
    '''
    GTDB-TK annotation of NC MAGs.
    '''

    gtdb_file = os.path.join(args.b,'nc_gtdbtk_MQNC/gtdbtk.bac120.summary.tsv')
    gtdb_file = '07_binannotation/bacteria/nc_gtdbtk_MQNC/taxonomy.txt''
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
            tax = line[1]
            
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
    fileout1 = os.path.join(args.o,'VAMBV3.MAGmotherbin.tax.txt')
    fileout1 = '10_abundances/VAMBV3.MAGmotherbin.tax.txt'

    with open(fileout1,'w') as out:
        out.write('motherbin' + '\t' + '\t'.join(lineage) + '\n')
        for MAG in MAG_tax:
            taxstring = [ MAG_tax[MAG][lin] for lin in lineage ]
            out.write(MAG + '\t' + '\t'.join(taxstring) + '\n')



if __name__ == "__main__":
    args = parser.parse_args()

    samples, sample_list, breadth_of_coverage, contig_order = get_samples_and_coverage(args)


    ### Subset original Martix to Viral Contigs
    bin_abundances, motherbins_avgRPKM, motherbins_avgRPM = get_bin_abundances(args,sample_list, breadth_of_coverage, contig_order)
    
    write_matrices(args,samples,  motherbins_avgRPKM, motherbins_avgRPM)
    
    read_write_MAG_taxonomy(args)