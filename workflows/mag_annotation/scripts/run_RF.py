#!/usr/bin/python
import sys
import argparse
import vambtools as _vambtools
import run_RF_modules
import collections as _collections
import os
import numpy as _np


parser = argparse.ArgumentParser(
    description="""Command-line benchmark utility.""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=True)

parser.add_argument('fastafile', help='Path to concatenated assembly VAMB')
parser.add_argument('clusterspath', help='Path to clusters.tsv')
parser.add_argument('annotationdir', help='Path to directory with contig annotations')
parser.add_argument('directoryout', help='Path to directory out')
parser.add_argument('-m', dest='min_bin_size', metavar='', type=int,
                    default=5000, help='Minimum size of bins - default [5000]')
parser.add_argument('-s', dest='separator', help='Binsplit separator', default=None)


def write_concat_bins(directory, bins, fastadict, compressed=False, maxbins=250, minsize=5000): 
    """Writes bins as FASTA files in a directory, one file per bin.
    Inputs:
        directory: Directory to create or put files in
        bins: {'name': {set of contignames}} dictionary (can be loaded from
        clusters.tsv using vamb.cluster.read_clusters)
        fastadict: {contigname: FastaEntry} dict as made by `loadfasta`
        compressed: Sequences in dict are compressed [False]
        maxbins: None or else raise an error if trying to make more bins than this [250]
        minsize: Minimum number of nucleotides in cluster to be output [0]
    Output: None
    """
    import os as _os
    import gzip as _gzip
    import vambtools as _vambtools

    import random

    # Safety measure so someone doesn't accidentally make 50000 tiny bins
    # If you do this on a compute cluster it can grind the entire cluster to
    # a halt and piss people off like you wouldn't believe.
    if maxbins is not None and len(bins) > maxbins:
        raise ValueError('{} bins exceed maxbins of {}'.format(len(bins), maxbins))

    # Check that the directory is not a non-directory file,
    # and that its parent directory indeed exists
    abspath = _os.path.abspath(directory)
    parentdir = _os.path.dirname(abspath)

    if parentdir != '' and not _os.path.isdir(parentdir):
        raise NotADirectoryError(parentdir)

    if _os.path.isfile(abspath):
        raise NotADirectoryError(abspath)

    if minsize < 0:
        raise ValueError("Minsize must be nonnegative")

    # Check that all contigs in all bins are in the fastadict
    allcontigs = set()

    for contigs in bins.values():
        allcontigs.update(set(contigs))

    allcontigs -= fastadict.keys()
    if allcontigs:
        nmissing = len(allcontigs)
        raise IndexError('{} contigs in bins missing from fastadict'.format(nmissing))

    # Make the directory if it does not exist - if it does, do nothing
    try:
        _os.mkdir(directory)
    except FileExistsError:
        pass
    
    bins_entries = []
    # Now actually print all the contigs to files
    for binname, contigs in bins.items():
        
        # Concatenate sequences of the bin
        concat_sequence = bytearray()
        for contig in contigs:
            entry = fastadict[contig]
            if compressed:
                uncompressed = bytearray(_gzip.decompress(entry.sequence))
                concat_sequence += uncompressed
            else:
                uncompressed = bytearray(entry.sequence)
                concat_sequence += uncompressed

        bin_entry = _vambtools.FastaEntry(binname, concat_sequence)           
        # Skip bin if it's too small
        if len(bin_entry.sequence) < minsize:
            continue
        bins_entries.append(bin_entry)
    
    random.shuffle(bins_entries)
    print('Writing:',len(bins_entries) ,'bins to file')
    filename = _os.path.join(directory, 'vamb_bins.1.fna')
    i = 1
    j = 1
    file = open(filename,'w')
    for entry in bins_entries:
        if i % 100000 == 0:
            j += 1 
            file.close()
            filename = _os.path.join(directory, 'vamb_bins.' + str(j) + '.fna')
            file = open(filename,'w')
        i += 1
        print(entry.format(),file=file)

def write_phamb_tables(RF_results,directory):
    '''Write Input table to RF-model and RF predictions'''
    try:
        os.mkdir(directory)
    except FileExistsError:
        pass

    annotation_table = os.path.join(directory,'vambbins_aggregated_annotation.txt')
    with open(annotation_table,'w') as out:
        header = ['binname','size','micomplete','VOG','dvf_score']
        out.write('\t'.join(header)+'\n')
        for i,binname in enumerate(RF_results.genome_order):
            row = [binname] + RF_results.df[i]
            out.write('\t'.join([str(i) for i in row])+'\n')

    prediction_table = os.path.join(directory,'vambbins_RF_predictions.txt')
    with open(prediction_table,'w') as out:
        header = ['binname','label','probability']
        out.write('\t'.join(header)+'\n')
        for i,binname in enumerate(RF_results.genome_order):
            row = RF_results.RF_predictions[i][1:]
            out.write('\t'.join([str(i) for i in row])+'\n')


### Relevant for the RF-model prediction
from sklearn.metrics import precision_score, recall_score, roc_auc_score, roc_curve, confusion_matrix, f1_score, auc, matthews_corrcoef
from scipy import sparse
import joblib

class RF_model():
    def __init__(self,RF_model,genomes):
        self.RF_model = RF_model
        self.RF_predictions = None
        self.RF_non_bacteria = None
        self.df = None
        self.genome_order = None
        
        ### Prepare and run PHAMB
        genome_order, sparse_df, regular_df = self._return_RF_prediction_dataframe(genomes)
        RF_predictions = self._run_RF_model(self.RF_model,genome_order,sparse_df)
        self.RF_predictions = RF_predictions
        self.genome_order = genome_order
        self.df = regular_df
        self.RF_non_bacteria = [row[1] for row in self.RF_predictions if row[2] == 'viral' ]

    @classmethod
    def _return_RF_prediction_dataframe(cls,genomes):
        '''On genome level: Return sparse-matrix for Prediction with PHAMB model'''
        annotation_types = ['micompletehmm','voghmm','deepvirfinder']
        df = []
        genome_order = []
        for genome in genomes.values():
            genome_order += [genome.name]
            row = [genome.totalsize]
            for type in annotation_types:
                value = 0
                if type in genome.genome_annotation:
                    if type == 'deepvirfinder':
                        value = genome.genome_annotation[type]['weighted_mean_score']
                    else:
                        value = genome.genome_annotation[type]
                row += [value]
            df.append(row)
        sparse_df = sparse.csr_matrix(df)
        return genome_order, sparse_df, df

    @classmethod
    def _run_RF_model(cls,RF_model,genome_order, sparse_df):
        '''ON GENOME LEVEL: Runs RF predictive-PHAMB model'''

        print('Loading Model and annotation table')
        trained_model = joblib.load(RF_model)

        predicted_genome_labels = trained_model.predict(sparse_df)
        prediction_probabilities = trained_model.predict_proba(sparse_df)
        predicted_genome_labels = [label.lower() for label in list(predicted_genome_labels) ] 
        rows = []
        for i, genome_name in enumerate(genome_order):
            rows.append(['PHAMB',genome_name ,predicted_genome_labels[i], prediction_probabilities[i][1]])
        return rows

    



if __name__=='__main__':

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    try:
        os.mkdir(args.directoryout)
    except FileExistsError:
        pass
    
    with open(args.clusterspath) as file:
        clusters = _vambtools.read_clusters(file)

    with _vambtools.Reader(args.fastafile, 'rb') as infile:
        fastadict = _vambtools.loadfasta(infile,compress=False)

    reference = run_RF_modules.Reference.from_clusters(clusters = clusters, fastadict=fastadict, minimum_contig_len=2000)
    annotation_directory = args.annotationdir
    viral_annotation_files = {
                    'deepvirfinder':os.path.join(annotation_directory,'all.DVF.predictions.txt'),
                    'voghmm':os.path.join(annotation_directory,'all.hmmVOG.tbl'),
                    'micompletehmm':os.path.join(annotation_directory,'all.hmmMiComplete105.tbl'),
                    }

    viral_annotation = run_RF_modules.Viral_annotation(annotation_files=viral_annotation_files,genomes=reference)

    rf_model_file = 'mag_annotation/dbs/RF_model.sav'
    RF_results = RF_model(rf_model_file,  genomes = viral_annotation.genomes)

    bins = {binname:clusters[binname] for binname in RF_results.RF_non_bacteria}
    write_concat_bins(os.path.join(args.directoryout,'vamb_bins'), bins, fastadict, compressed=False, maxbins=len(bins), minsize=args.min_bin_size)
    write_phamb_tables(RF_results,args.directoryout)