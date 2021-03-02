from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys 
import os

class bin_annotation:
    def __init__(self,bin_name,motherbin):
        self.motherbin = motherbin
        self.bin_name = bin_name
        self.ncontigs = None
        self.binsize = 0


def read_in_clusters(cluster_file):
    
    print('Reading in Clusters, this may take a while...')
    bin_annotations = dict()
    clusters = dict()
    with open(cluster_file,'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            cluster, contig = line[0],line[1]
            sample = contig.split('_')[0]
            binid = sample + '_' + cluster
            if 'length' in contig:
                contiglength = int(contig.split('length')[1].split('_')[1])
            if binid not in bin_annotations:
                bin_annotations[binid] = bin_annotation(binid,cluster)
                bin_annotations[binid].ncontigs = 1
                bin_annotations[binid].binsize += contiglength
            else:
                bin_annotations[binid].ncontigs += 1
                bin_annotations[binid].binsize += contiglength
            if contig not in clusters:
                clusters[contig] = (cluster,binid)
    return clusters, bin_annotations


def parse_MVX_hits(MVX_file,clusters,bin_annotations,seqid=95,coverage = 0.5, NC_MVX=False ):
    '''
    Calculate sequence coverage of each MVX genome/sequence by each MGX bin 
    Also keep track of how many contigs in a MGX bin maps to the given MVX sequence.
    '''
    MVX_genomes_coverage = dict()
    bin_MVXhits = dict()
    with open(MVX_file,'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            contig = line[0]
            sequence_identity = float(line[2])
            alignment_length = int(line[3])
            MVX_genome = line[1]
            MVX_genome_length = int(line[-1])
            if contig in clusters:
                binid = clusters[contig][1]
                if sequence_identity >= seqid and alignment_length/MVX_genome_length >= coverage:     
                    if binid not in bin_MVXhits:
                        bin_MVXhits[binid] = set([contig])
                    else:
                        bin_MVXhits[binid].add(contig)
                        
                    if MVX_genome not in MVX_genomes_coverage:
                        MVX_genomes_coverage[MVX_genome] = np.zeros(MVX_genome_length)
                    
                    ref_from = int(line[7])
                    ref_to = int(line[8])
                    if ref_to > ref_from:
                        MVX_genomes_coverage[MVX_genome][ref_from:ref_to] += 1
                    else:
                        MVX_genomes_coverage[MVX_genome][ref_to:ref_from] += 1
    
    ### Calculate MVX genome coverage
    for MVX_genome in MVX_genomes_coverage:
        coverage = MVX_genomes_coverage[MVX_genome]
        covered = 1-(len(coverage[coverage==0])/len(coverage))
        MVX_genomes_coverage[MVX_genome] = covered
    MVX_genome_df = pd.DataFrame.from_dict(MVX_genomes_coverage, orient='index',columns=['coverage'])
    completeness_fact = {'LQ':25, 'MQ':50, 'NC':90}
    completeness = []
    for coverage in MVX_genome_df['coverage']:
        percentage_coverage = coverage*100
        decision = 'NA'
        for quality in completeness_fact:
            if percentage_coverage >= completeness_fact[quality]:
                decision = quality
        completeness.append(decision)
    MVX_genome_df['Q'] = completeness
    
    ### Summarise strong hits for each bin
    bin_MVX_fraction = dict()
    for binid in bin_MVXhits:
        n_contigs = bin_annotations[binid].ncontigs
        nhits = len(bin_MVXhits[binid])
        bin_imgvr_fraction = nhits/n_contigs
        bin_MVX_fraction[binid] = bin_imgvr_fraction

    return bin_MVX_fraction, MVX_genome_df


def prepare_train_test_set(vamb_bins_file,clusters_file,MVX_blast_file,checkm_file):
    '''
    Each VAMB bin is annotated either as Bacterial or Viral.
    - Bacteria if CheckM completeness >= 25% and contamination <= 10%
    - Viral if 95% of all contigs in MGX bin maps to an MVX genome with minimum 95% sequence identity and minimum 50% query sequence coverage.
    '''
    vambbins_df = pd.read_csv(vamb_bins_file,sep='\t')
    LABEL = ['NA' for i in range(vambbins_df.shape[0])]

    ### Assign Bacterial Labels using CheckM result
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
            if float(completeness) >= 10 and float(contamination) <= 30:
                MQNC_bins.add(binid)

    ### Bacterial first
    for i,binid in enumerate(vambbins_df['binid']):
        if binid in MQNC_bins:
            LABEL[i] = 'Bacterial'
    
    ### Viral     
    clusters, binsannos = read_in_clusters(cluster_file = clusters_file)
    bin_MVX_fraction,MVX_genome_df = parse_MVX_hits(MVX_blast_file, clusters, binsannos)    
    nhallm_dict = dict(zip(vambbins_df.binid, vambbins_df.nhallm))
    VOG_dict = dict(zip(vambbins_df.binid, vambbins_df.nVOGs))
    DVF_dict = dict(zip(vambbins_df.binid, vambbins_df.cluster_DVF_score))
    for i,binid in enumerate(vambbins_df['binid']):
        if binid in bin_MVX_fraction:
            mvx_score = bin_MVX_fraction[binid]
            if mvx_score >= 0.95:
                LABEL[i] = 'Viral'
            
    ### Remove all bins without an Annotation
    vambbins_df['LABEL'] = LABEL
    vambbins_df_subset = vambbins_df[['binid','binsize','nhallm','nVOGs','cluster_DVF_score','LABEL']]
    bacterial_mask = list( (vambbins_df_subset.LABEL=='NA') & (vambbins_df_subset.nhallm>95) )
    vambbins_df_subset.loc[bacterial_mask,'LABEL'] = 'Bacterial'
    vambbins_df_subset = vambbins_df_subset[vambbins_df_subset.LABEL!='NA']
    labels = list(vambbins_df_subset.LABEL)

    # Split 60% of annotated bins into test and 40% for the test
    train, test, train_labels, test_labels = train_test_split(vambbins_df_subset, labels, 
                                                          stratify = labels,
                                                          test_size = 0.6, 
                                                          random_state = 1)
    train_objects = {'train':train, 'test':test, 'train_labels':train_labels, 'test_labels':test_labels}
    return vambbins_df, vambbins_df_subset ,MVX_genome_df, train_objects




def randomsearch_crossval(training_set, training_labels):
    """[summary]

    Args:
        training_set ([type]): [description]
        training_labels ([type]): [description]
    """

    from sklearn.ensemble import RandomForestClassifier
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import RandomizedSearchCV
    from sklearn.model_selection import GridSearchCV
    rf = RandomForestRegressor(random_state = 42)
    from pprint import pprint

    # Number of trees in random forest
    n_estimators = [int(x) for x in np.linspace(start = 10, stop = 500, num = 20)]
    # Number of features to consider at every split
    max_features = ['auto', 'sqrt']
    # Maximum number of levels in tree
    max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
    max_depth.append(None)
    # Minimum number of samples required to split a node
    min_samples_split = [2, 5, 10]
    # Minimum number of samples required at each leaf node
    min_samples_leaf = [1, 2, 4]
    # Method of selecting samples for training each tree
    bootstrap = [True, False]
    # Create the random grid
    random_grid = {'n_estimators': n_estimators,
                'max_features': max_features,
                'max_depth': max_depth,
                'min_samples_split': min_samples_split,
                'min_samples_leaf': min_samples_leaf,
                'bootstrap': bootstrap}

    default_model = RandomForestClassifier(n_estimators=300, 
                            random_state=1, 
                            max_features = 'sqrt',
                            n_jobs=1, verbose = 1)

    rf_random = RandomizedSearchCV(estimator = default_model, param_distributions = random_grid, 
                               n_iter = 50, cv = 3, verbose=2, random_state=42, n_jobs = 1)

    # Fit the random search model
    rf_random.fit(np.array(training_set), np.array(training_labels))
    best_random_model = rf_random.best_estimator_


    return(best_random_model)

def gridsearch_crossval(training_set, training_labels):
    """[summary]

    Returns:
        [type]: [description]
    """


    from sklearn.ensemble import RandomForestClassifier
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import RandomizedSearchCV
    from sklearn.model_selection import GridSearchCV
    rf = RandomForestRegressor(random_state = 42)
    from pprint import pprint


    param_grid = {
    'bootstrap': [True],
    'max_depth': [20, 30,40],
    'max_features': ['auto'],
    'min_samples_leaf': [1,2],
    'min_samples_split': [2,4],
    'n_estimators': [100,345,400]
    }
    default_model = RandomForestClassifier(n_estimators=300, 
                            random_state=1, 
                            max_features = 'sqrt',
                            n_jobs=1, verbose = 1)


    grid_search = GridSearchCV(estimator = default_model, param_grid = param_grid, 
                          cv = 3, n_jobs = 1, verbose = 2)
    # Fit the grid search to the data
    grid_search.fit(np.array(training_set), np.array(training_labels))

    best_grid_model = grid_search.best_estimator_
    return(best_grid_model)




### Code to return contig table to assert DVF performance 
def return_DVF_table(cluster_file, DVF_predictions_file, labelled_bins):
    '''Function to parse DVF predictions for single contigs
       This is used for assessing performance of Viral prediction if it was conducted without Bins
    '''
    clusters, binsannos = read_in_clusters(cluster_file)
    DVF_table = []
    
    ### Parse DVF predictions
    with open(DVF_predictions_file,'r') as infile:
        for line in infile:
            contig, score,pvalue, sample = line.strip().split('\t')
            if contig not in clusters:
                continue
            binid = clusters[contig][1]
            if binid not in labelled_bins:
                continue 
            true_label = labelled_bins[binid]
            prediction = None

            ### DVF prediction determined a contig as Viral if both P-value < 0.05 and Score above 0.5 
            if float(pvalue) < 0.05 and float(score) > 0.5:
                prediction = 'Viral'
            else:
                prediction = 'Bacterial'
            row = [contig,binid,true_label,prediction,score]
            DVF_table.append(row)
    DVF_table = pd.DataFrame(DVF_table,columns=['contigid','binid','LABEL','DVF_prediction','DVF_score'])
    return DVF_table



### Functions for plotting
def train_and_test_RF(model,train, test, train_labels, test_labels):

    model.fit(train, train_labels)

    train_predictions = model.predict(train)
    train_probs = model.predict_proba(train)[:, 1]
    print(f'Train ROC AUC Score: {roc_auc_score(train_labels, train_probs)}')


    ### Test Set
    predictions = model.predict(test)
    probs = model.predict_proba(test)[:, 1]

    print(f'Test ROC AUC Score: {roc_auc_score(predictions, probs)}')
    return train_predictions,train_probs, predictions, probs, model

def plot_precision_recall(predictions,probs,bin_test_labels):
    
    lr_precision, lr_recall, _ = precision_recall_curve(bin_test_labels, probs)
    bin_predictions = np.where(predictions == 'Bacterial',0,1)
    lr_f1, lr_auc = f1_score(bin_test_labels, bin_predictions), auc(lr_recall, lr_precision)
        
    
    # summarize scores
    print('Logistic: f1=%.3f auc=%.3f' % (lr_f1, lr_auc))
    # plot the precision-recall curves
    no_skill = len(bin_test_labels[bin_test_labels==1]) / len(bin_test_labels)
    plt.plot([0, 1], [no_skill, no_skill], linestyle='--', label='No Skill')
    plt.plot(lr_recall, lr_precision, marker='.', label='RF-Model')
    # axis labels
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    # show the legend
    plt.legend()
    # show the plot
    plt.show()


def plot_precision_recall_wDVF(predictions,probs,bin_test_labels,DVF_probs,DVF_labels):
    
    lr_precision, lr_recall, _ = precision_recall_curve(bin_test_labels, probs)
    bin_predictions = np.where(predictions == 'Bacterial',0,1)
    lr_f1, lr_auc = f1_score(bin_test_labels, bin_predictions), auc(lr_recall, lr_precision)
    
    ### DVF 
    
    DVF_lr_precision, DVF_lr_recall, _ = precision_recall_curve(DVF_labels, DVF_probs)
    DVF_predictions = np.where(DVF_labels == 'Bacterial',0,1)
    
    
    
    # summarize scores
    print('Logistic: f1=%.3f auc=%.3f' % (lr_f1, lr_auc))
    # plot the precision-recall curves
    no_skill = len(bin_test_labels[bin_test_labels==1]) / len(bin_test_labels)
    plt.plot([0, 1], [no_skill, no_skill], linestyle='--', label='No Skill')
    plt.plot(lr_recall, lr_precision, marker='.', label='RF-Model')
    plt.plot(DVF_lr_recall, DVF_lr_precision, marker='.', label='DVF-bin-score')
    # axis labels
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    # show the legend
    plt.legend()
    # show the plot
    plt.show()

def plot_AUC_curve(probs,bin_test_labels):
    
    ### Make Random Guess array
    ns_probs = [0 for _ in range(len(bin_test_labels))]
    # calculate scores
    ns_auc = roc_auc_score(bin_test_labels, ns_probs)
    lr_auc = roc_auc_score(bin_test_labels, probs)
    # summarize scores
    print('No Skill: ROC AUC=%.3f' % (ns_auc))
    print('RF-model: ROC AUC=%.3f' % (lr_auc))
    # calculate roc curves
    ns_fpr, ns_tpr, _ = roc_curve(bin_test_labels, ns_probs)
    lr_fpr, lr_tpr, _ = roc_curve(bin_test_labels, probs)
    # plot the roc curve for the model
    plt.plot(ns_fpr, ns_tpr, linestyle='--', label='No Skill')
    plt.plot(lr_fpr, lr_tpr, marker='.', label='RF-Model')
    # axis labels
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # show the legend
    plt.legend()
    # show the plot
    plt.show()
    

def plot_AUC_curve_wDVF(probs,bin_test_labels, DVF_probs,DVF_labels):
    
    ### Make Random Guess array
    ns_probs = [0 for _ in range(len(bin_test_labels))]
    # calculate scores
    ns_auc = roc_auc_score(bin_test_labels, ns_probs)
    lr_auc = roc_auc_score(bin_test_labels, probs)
    # summarize scores
    print('No Skill: ROC AUC=%.3f' % (ns_auc))
    print('RF-model: ROC AUC=%.3f' % (lr_auc))
    # calculate roc curves
    ns_fpr, ns_tpr, _ = roc_curve(bin_test_labels, ns_probs)
    lr_fpr, lr_tpr, _ = roc_curve(bin_test_labels, probs)
    
    DVF_lr_fpr, DVF_lr_tpr, _ =  roc_curve(DVF_labels, DVF_probs)
    # plot the roc curve for the model
    plt.plot(ns_fpr, ns_tpr, linestyle='--', label='No Skill')
    plt.plot(lr_fpr, lr_tpr, marker='.', label='RF-Model')
    plt.plot(DVF_lr_fpr, DVF_lr_tpr, marker='*', label='DVF-bin-score')
    # axis labels
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # show the legend
    plt.legend()
    # show the plot
    plt.show()
