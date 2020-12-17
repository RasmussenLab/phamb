#!/bin/Python
### Python code to parse .npz into depth matrix in .txt format with contignames and contiglengths
import numpy as np
import sys 
import argparse
import os
import scipy.stats
import multiprocessing as mp
import pickle

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
   ''')
parser.add_argument('-d', help='JGI depth matrix file')
parser.add_argument('-l', help='List of JGI coverage files')
parser.add_argument('-m', help='checkm-mags file')
parser.add_argument('-c', help='Cluster file - will be used to extract specific contigs of interest',default=None)
parser.add_argument('-o', help='file outbase')



def get_samples_and_coverage(args):
    jgi_samples = args.l
    jgi_samples = '05_binning/vamb_on_jgi_v3/HMP2/HMP2.depthfiles.complete.txt'
    sample_list = []
    sample_counter = 0
    with open(jgi_samples,'r') as infile:
        samples = ['contigname']
        for line in infile:
            line = line.strip()
            base = os.path.basename(line)
            base = base.replace('.jgi.depth.txt','')
            coverage_file = line.replace('.jgi.depth.txt','.position.zero.coverage')
            if os.path.exists(coverage_file):
                sample_list.append((coverage_file,base,sample_counter))
                samples.append(base)
            sample_counter += 1
    samples = np.array(samples).reshape(1,-1)

    ### Get the right contig order from any given JGI file
    jgi_depth_file = sample_list[0][0].replace('.position.zero.coverage','.jgi.depth.txt')
    contig_order = []
    with open(jgi_depth_file,'r') as infile:
        infile.readline()
        for line in infile:
            contig = line.strip().split()[0]
            contig_order.append(contig)

    ### Calculate breadth of coverage for each contig in each sample 
    outbase = '10_abundances'
    if not os.path.exists(outbase + '/contig_breadth_of_coverage.pickle'):

        def get_coverage_worker(coverage_file,sample,return_dict):
            '''
            For each contig, contig coverage is determined by the proportion of non-covered positions. 
            Notice we only consider how large a proportion the '0' depth position accounts for.
            '''

            breadth_of_coverage = dict()
            with open(coverage_file,'r') as infile:
                for line in infile:
                    contig, depth, startpos, endpos, fraction = line.strip().split('\t')
                    if depth == '0':
                        contig_coverage = 1-float(fraction)
                        if contig_coverage >= 0.75:
                            breadth_of_coverage[contig] = 1
            return_dict[sample] = breadth_of_coverage
            
        ### Using Multiple threads!!! To parse the coverage files for each contig.
        manager = mp.Manager()
        return_dict =  manager.dict()
        n_cores = 3
        pool = mp.Pool(n_cores)
        for item in sample_list:
            coverage_file = item[0]
            sample = item[1]
            pool.apply_async(get_coverage_worker, args=(coverage_file,sample,return_dict))
        pool.close()
        pool.join()
        breadth_of_coverage = dict() 
        for sample in return_dict:
            breadth_of_coverage[sample] = return_dict[sample]
        
        ### Save the breadth of coverage
        outbase = '10_abundances'
        with open(outbase + '/contig_breadth_of_coverage.pickle', 'wb') as fp:
            pickle.dump(breadth_of_coverage, fp, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        breadth_of_coverage = pickle.load( open( outbase + '/contig_breadth_of_coverage.pickle', "rb" ) )
      
    return samples, sample_list, breadth_of_coverage, contig_order



class bin_abundance:
     def __init__(self,bin_name,mothercluster):
        self.bin = bin_name
        self.binsize = 0
        self.mothercluster = mothercluster
        self.contigs = set()
        self.binRPKM = None
        self.binRPM = None
     


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
        for line in infile:
            line = line.strip().split('\t')[0]
            mags.add(line)

    bin_abundances = dict()
    with open(clusterfile,'r') as infile:
        for line in infile:
            line = line.strip().split()
            mothercluster, contig = line[0],line[1]
            binid = contig.split('_')[0] + '_' + mothercluster
            if not binid in mags:
                continue
            contiglen = contig_lengths_lookup[contig]
            if binid not in bin_abundances:     
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