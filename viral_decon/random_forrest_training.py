#!/bin/python

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import os
import pandas as pd 
import numpy as np
from sklearn.metrics import precision_score, recall_score, roc_auc_score, roc_curve, confusion_matrix



### https://github.com/WillKoehrsen/Machine-Learning-Projects/blob/master/Random%20Forest%20Tutorial.ipynb




### Viral Next
### Only for Training and teseting


class bin_annotation:
    def __init__(self,bin_name,motherbin):
        self.motherbin = motherbin
        self.bin_name = bin_name
        self.ncontigs = None
        self.binsize = 0

def read_in_clusters(cluster_file):
    
    print('Reading in Clusters, this may take a while...')
    binsannos = dict()
    cls = dict()
    with open(cluster_file,'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            cluster, contig = line[0],line[1]
            sample = contig.split('_')[0]
            binid = sample + '_' + cluster
            if 'length' in contig:
                contiglength = int(contig.split('length')[1].split('_')[1])


            if binid not in binsannos:
                binsannos[binid] = bin_annotation(binid,cluster)
                binsannos[binid].ncontigs = 1
                binsannos[binid].binsize += contiglength

            else:
                binsannos[binid].ncontigs += 1
                binsannos[binid].binsize += contiglength

            if contig not in cls:
                cls[contig] = (cluster,binid)
    return cls, binsannos



def parse_MVX_hits(MVX_file,cls,binsannos,seqid=95,coverage = 0.5):
    '''
    Calculate How many Contigs in each bin with hits to IMGVR genomes
    '''

    bin_MVXhits = dict()
    with open(MVX_file,'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            contig = line[0]
            sequence_identity = float(line[2])
            alignment_length = int(line[3])
            imgvr_genome_length = int(line[-1])
            if contig in cls:
                binid = cls[contig][1]

                if sequence_identity >= seqid and alignment_length/imgvr_genome_length >= coverage:
                    if binid not in bin_MVXhits:
                        bin_MVXhits[binid] = set([contig])
                    else:
                        bin_MVXhits[binid].add(contig)
    
    ### Summarise strong hits for each bin
    bin_MVX_fraction = dict()
    for binid in bin_MVXhits:
        n_contigs = binsannos[binid].ncontigs
        nhits = len(bin_MVXhits[binid])
        bin_imgvr_fraction = nhits/n_contigs
        bin_MVX_fraction[binid] = bin_imgvr_fraction

    return bin_MVX_fraction


def get_MVX_bins(cluster_file,MVX_file):
    cls, binsannos = read_in_clusters(cluster_file = cluster_file)
    bin_MVX_fraction = parse_MVX_hits(MVX_file, cls, binsannos)
    return bin_MVX_fraction


def prepare_train_test_set(vamb_bins_file,clusters_file, MVX_file,checkm_file):
    vambbins_df = pd.read_csv(vamb_bins_file,sep='\t')

    ### Assign Bacterial Labels
    LABEL = ['NA' for i in range(vambbins_df.shape[0])]

    # Using CHECKM

    MQNC_bins = set()
    with open(checkm_file,'r') as infile:
        for line in infile:
            if line[0] =='-':
                continue
            elif 'Bin Id' in line:
                continue
            line = line.strip().split()
            line = [line[0]]+[' '.join(line[1:3])] + line[3:]
            binid = line[0]
            completeness = float(line[-3])
            contamination = float(line[-2])
            if float(completeness) >= 25 and float(contamination) <= 10:
                MQNC_bins.add(binid)

    ### Bacterial firsst
    for i,binid in enumerate(vambbins_df['binid']):
        if binid in MQNC_bins:
            LABEL[i] = 'Bacterial'
    
    ### Viral 
    bin_MVX_fraction = get_MVX_bins(clusters_file,MVX_file)

    for i,binid in enumerate(vambbins_df['binid']):
        if binid in bin_MVX_fraction:
            mvx_score = bin_MVX_fraction[binid]
            if mvx_score >= 0.5:
                LABEL[i] = 'Viral'

    ### Remove all bins without an Annotation
    vambbins_df['LABEL'] = LABEL

    return vambbins_df




def evaluate_model(predictions, probs, train_predictions, train_probs):
    """Compare machine learning model to baseline performance.
    Computes statistics and shows ROC curve.
    """
    
    baseline = {}
    
    baseline['recall'] = recall_score(test_labels, [1 for _ in range(len(test_labels))])
    baseline['precision'] = precision_score(test_labels, [1 for _ in range(len(test_labels))])
    baseline['roc'] = 0.5
    
    results = {}
    
    results['recall'] = recall_score(test_labels, predictions)
    results['precision'] = precision_score(test_labels, predictions)
    results['roc'] = roc_auc_score(test_labels, probs)
    
    train_results = {}
    train_results['recall'] = recall_score(train_labels, train_predictions)
    train_results['precision'] = precision_score(train_labels, train_predictions)
    train_results['roc'] = roc_auc_score(train_labels, train_probs)
    
    for metric in ['recall', 'precision', 'roc']:
        print(f'{metric.capitalize()} Baseline: {round(baseline[metric], 2)} Test: {round(results[metric], 2)} Train: {round(train_results[metric], 2)}')


def train_and_test_RF(model,train, test, train_labels, test_labels):

    model.fit(train, train_labels)

    train_predictions = model.predict(train)
    train_probs = model.predict_proba(train)[:, 1]
    print(f'Train ROC AUC Score: {roc_auc_score(train_labels, train_probs)}')


    ### Test Set
    predictions = model.predict(test)
    probs = model.predict_proba(test)[:, 1]

    print(f'Test ROC AUC Score: {roc_auc_score(predictions, probs)}')


    results = {}
    results['recall'] = recall_score(test_labels, predictions,average='macro')
    results['precision'] = precision_score(test_labels, predictions,average='macro')
    results['roc'] = roc_auc_score(test_labels, probs,average='macro')

    train_results = {}
    train_results['recall'] = recall_score(train_labels, train_predictions,average='macro')
    train_results['precision'] = precision_score(train_labels, train_predictions,average='macro')
    train_results['roc'] = roc_auc_score(train_labels, train_probs,average='macro')
    return train_results, results





########## Random Forest CODE and stuff ########## 

### Diabimmune 
vamb_bins_file = '../diabimu_t1d/04_annotation/annotation_summaries_new/vambbins_aggregated_annotation.txt'
clusters_file = '../diabimu_t1d/05_binning/vamb_on_jgi_v3/Diabimmune/clusters.tsv'
MGX_MVX_blast = '../diabimu_t1d/06_blast_MGXMVX/MGX.2000.combined.diabimmune_viral_contigs.m6'
checkm_file ='../diabimu_t1d/07_binannotation/bacteria/checkm/Diabimmune.output'
diab_vambbins_df = prepare_train_test_set(vamb_bins_file,clusters_file, MGX_MVX_blast,checkm_file)
diab_subset = diab_vambbins_df[ diab_vambbins_df['LABEL'] != 'NA' ] 
evaluation_labels = np.array(diab_subset.pop('LABEL'))
diab_subset =  diab_subset[['binsize','nhallm','nVOGs','cluster_DVF_score']]



### COPSAC
vamb_bins_file = '04_annotation/annotation_summaries/vambbins_aggregated_annotation.txt'
clusters_file = '05_binning/vamb_on_jgi/COPSAC/clusters.tsv'
MGX_MVX_blast = '06_blast_MGXMVX/MGX.2000.copsac.assemblies.vOTUs.selected.m6'
checkm_file ='07_binannotation/bacteria/checkm/COPSAC.output'
vambbins_df = prepare_train_test_set(vamb_bins_file,clusters_file, MGX_MVX_blast, checkm_file)
subset = vambbins_df[ vambbins_df['LABEL'] != 'NA' ] 
 
labels = np.array(subset.pop('LABEL'))
subset = subset[['binsize','nhallm','nVOGs','cluster_DVF_score']]


# 30% examples in test data
train, test, train_labels, test_labels = train_test_split(subset, labels, 
                                                          stratify = labels,
                                                          test_size = 0.3, 
                                                          random_state = 1)
model = RandomForestClassifier(n_estimators=100, 
                            random_state=1, 
                            max_features = 'sqrt',
                            n_jobs=1, verbose = 1)

### Train on COPSAC
train_results, test_results = train_and_test_RF(model, train, test, train_labels, test_labels )


### Evaluate on Diabimmune Test-set
eval_predictions = model.predict(diab_subset)
eval_probs = model.predict_proba(diab_subset)[:, 1]

print(f'Evaluation ROC AUC Score: {roc_auc_score(evaluation_labels, eval_probs)}')
eval_results = {}
eval_results['recall'] = recall_score(evaluation_labels, eval_predictions,average='macro')
eval_results['precision'] = precision_score(evaluation_labels, eval_predictions,average='macro')
eval_results['roc'] = roc_auc_score(evaluation_labels, eval_probs,average='macro')

### 
confusion_matrix(evaluation_labels,eval_predictions)

