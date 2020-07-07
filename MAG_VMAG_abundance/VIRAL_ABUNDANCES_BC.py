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
parser.add_argument('-l', help='List with paths to JGI coverage files')
parser.add_argument('-f', help='fasta index file')
parser.add_argument('-c', help='Cluster file - will be used to extract specific contigs of interest')
parser.add_argument('-o', help='file outbase')



def get_samples_and_coverage(args):
    jgi_samples = args.l
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
                        breadth_of_coverage[contig] = contig_coverage
        return_dict[sample] = breadth_of_coverage
        
    ### Using Multiple threads!!! To parse the coverage files for each contig.
    manager = mp.Manager()
    return_dict =  manager.dict()
    n_cores = 12
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
      
    return samples, sample_list, breadth_of_coverage, contig_order



class bin_abundance:
     def __init__(self,bin_name,mothercluster):
        self.bin = bin_name
        self.mothercluster = mothercluster
        self.contigs = set()
        self.binabundance = None
     


def get_bin_abundances(args,sample_list, breadth_of_coverage, contig_order):

    ### parse contigs of cluster file 
    clusterfile = args.c
    
    bin_abundances = dict()
    viral_contigs = []
    clusterfile = 'analysis/vContact2bins/VAMBV3.Viral_RF_predictions.contigs.txt'
    with open(clusterfile,'r') as infile:
        for line in infile:
            line = line.strip().split()
            binid, contig = line[0],line[1]
            viral_contigs.append(contig)
            if binid not in bin_abundances:
                mothercluster = binid.split('_')[1]
                binobject = bin_abundance(binid,mothercluster)
                binobject.contigs.add(contig)
                bin_abundances[binid] = binobject
            else:
                binobject = bin_abundances[binid]
                binobject.contigs.add(contig)
    
    ### Parse depth matrix and calculate weighted mean
    matrixfile = args.d
    matrix = np.load(matrixfile)['arr_0']

    matrix_indices_dict = {}
    for i in range(len(contig_order)):
        contigname = contig_order[i]
        matrix_indices_dict[contigname] = i
    sample_indices = np.array([x[2] for x in sample_list])
    sample_order = [x[1] for x in sample_list]
    nsamples = len(sample_order)
    matrix_subset = matrix[:,sample_indices]

    ### Calculate Bin-abundance across all samples - if coverage of a contig is too low, set depth to zero.
    for binid in bin_abundances:
        binobject = bin_abundances[binid]
        bin_contigs = binobject.contigs

        contig_indices = [ matrix_indices_dict[contigname] for contigname in bin_contigs]
        binabundances = matrix_subset[np.array(contig_indices),]
        
        ### Make sure to set the JGI depth to zero in samples with low coverage for each contig, 
        # If contig in the dict breadth_of_coverage, its > 75% set to False.
        contig_indices = []
        for j,contig in enumerate(bin_contigs):
            ones_boolean = np.ones((1, nsamples) , dtype = bool)[0]
            for i,sample in enumerate(sample_order):
                if contig in breadth_of_coverage[sample]:
                    ones_boolean[i] = False
            binabundances[j,ones_boolean] = 0

        #bin_weighted_sample_averages = np.average(binabundances,axis=0, weights=contig_lengths)
        bin_sample_sums = np.sum(binabundances,axis=0)
        binobject.binabundance = bin_sample_sums
        bin_abundances[binid] = binobject
    
    ### Organise bins belonging to each Mothercluster
    motherbins = dict()
    for binid in bin_abundances:
        clustername = bin_abundances[binid].mothercluster
        if not clustername in motherbins:
            motherbins[clustername] = set([binid])
        else:
            motherbins[clustername].add(binid)
    
    ### Calculate a combined mothercluster abundance
    motherbins_abundances = {}
    for clustername in motherbins:
        matrixsubset = []
        bins = motherbins[clustername]
        for binid in bins:
            if binid in bin_abundances:
                matrixsubset.append(bin_abundances[binid].binabundance)
        
        ### Sum binabundances 
        clusterabundances = np.sum(matrixsubset, axis=0)
        motherbins_abundances[clustername] = clusterabundances


    ### Remember bin_abundances is a dict-object of class-objects
    ### motherbins_abundances is a standard dict 

    return bin_abundances, motherbins_abundances



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




def write_matrices(args ,samples ,bin_abundances, motherbins_abundances):
    
    outbase = args.o
    
    ### Write out Bin-abundances - heavy duty file...
    #fileout1 = os.path.join(outbase,'VAMBV3.bin.breadthcov.abundances.txt')
    #matrix = []
    #bins = []
    #for binid in bin_abundances:
    #    bins.append(binid)
    #    matrix.append(bin_abundances[binid].binabundance)
    #matrix = np.array(matrix)
    #matrix = matrix.astype('float32').round(2)
    #bins = np.array(bins).reshape(-1,1)
    #matrix = np.column_stack((bins,matrix))
    #header = ['bin'] + list(samples[0][1:])
    #header = '\t'.join(header)
    #np.savetxt(fileout1, matrix , delimiter='\t', header=header, newline='\n', fmt='%s')



    ### Write out Mothercluster Abundances
    fileout1 = os.path.join(outbase,'VAMBV3.motherbin.breadthcov.abundances.txt')
    mbin_median = [ (cluster,np.median(motherbins_abundances[cluster])) for cluster in motherbins_abundances]
    mbin_median = sorted(mbin_median , key= lambda x: x[1],reverse=True)


    matrix = []
    clusters = []
    for tup in mbin_median:
        cluster, val = tup
        clusters.append(cluster)
        matrix.append(motherbins_abundances[cluster])
    matrix = np.array(matrix)
    matrix = matrix.astype('float32').round(2)

    clusters = np.array(clusters).reshape(-1,1)
    matrix = np.column_stack((clusters,matrix))
    header = ['mothercluster'] + list(samples[0][1:])
    header = '\t'.join(header)
    np.savetxt(fileout1, matrix , delimiter='\t', header=header, newline='\n', fmt='%s')




if __name__ == "__main__":
    args = parser.parse_args()

    samples, sample_list, breadth_of_coverage, contig_order = get_samples_and_coverage(args)

    ### Subset original Martix to Viral Contigs
    bin_abundances, motherbins_abundances = get_bin_abundances(args,sample_list, breadth_of_coverage, contig_order)

    
    write_matrices(args,samples, bin_abundances, motherbins_abundances)