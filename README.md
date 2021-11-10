# phamb
A Phage from metagenomic bins (phamb) discovery approach used to isolate metagenome derived viromes and High-quality viral genomes.

The repository contains scripts and workflows used in our viral follow up study on the binning tool [VAMB](https://github.com/RasmussenLab/vamb) where we have benchmarked not only the quality and quantity of viral MAGs but also the viral overlap with metaviromes. In our analysis, [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) has been important for assessing the actual gain of using viral MAGs relative to single-contig evaluation, a big kudos to Nayfach et al. for this great tool. 

We have applied this approach to 3 different datasets and recovered up to 6,077 High-quality genomes from 1,024 viral populations, this is 200% more compared to only evaluation single-contigs. Similar to what we have observed for Bacterial bins, VAMB achieves high intra-VAMB-cluster ANI (>97.5%) also for viral bins, our best example here is accurate clustering of crAss-like bins found in the IBD Human Microbiome Project 2 dataset. 

- Our (recommended) workflow is to isolate the virome search space prior to running viral evaluation/prediction tools. For this, we have trained a Random Forest model on viral bins established using paired metagenomic and metavirome datasets. This massively helps in reducing computational time especially on larger datasets.
- We strongly advise to only use Medium-quality and High-quality viral bins evaluated using the AAI-model in CheckV. We found the HMM-model is not currently well-suited for viral MAGs. Low-quality viral bins may likely represent fragmented/incomplete viruses or general contamination, we advice caution with using these. 

## [Prerequisites] - Snakemake 

In order to run most of it, you need conda installed in order to install snakemake. Most snakemake workflows comes with conda-environments, thus dependencies and programmes are automatically installed. 

```bash
conda install -n snakemake_env snakemake pygraphviz python=3.8 cython scikit-learn==0.21.3

```

## 1. MAG annotation for isolating Metagenomic derived viromes

### Database and file requirements
VAMB bins and concatenated assemblies. 

```bash
contigs.fna.gz #Concatenated assembly 
vamb/clusters.tsv   #Clustered contigs based on the above contigs.fna.gz file 
```

Furthermore. 
* [VOGdb](https://vogdb.csb.univie.ac.at/download) - untar `vog.hmm.tar.gz` to get all hmm files. Concatenate all the hmm-files these into an `AllVOG.hmm` file  [File path needs to be specified in `config.yaml`]
* [Micomplete Bacterial HMMs](https://bitbucket.org/evolegiolab/micomplete/src/master/micomplete/share/Bact105.hmm)   [File path needs to be specified in `config.yaml`]
* Clone [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) git clone https://github.com/jessieren/DeepVirFinder


### How to Run 

Copy the phamb repository, extract the `mag_annotation` workflow and split contigs (using the provided script) to allow annotation to be run in parallel.

```bash
mkdir -p projectdir 
cd projectdir 
git clone the repository https://github.com/RasmussenLab/phamb.git
cp -r phamb/workflows/mag_annotation .
python mag_annotation/scripts/split_contigs.py -c contigs.fna.gz 

```

- Now the `contigs.fna.gz` is splitted into individual assemblies i.e. `assembly/{sample}/{sample}.fna`
- In addition, a `sample_table.txt` file is created with a line for each sample. Check that `sample_table.txt` contains sample identifiers corresponding to the ones you expect. The number of lines should correspond to the number of samples used to make the concatenated assembly (`contigs.fna.gz`).
- Now, specify paths for databases, vamb directory, location of assembly and computational resources in `mag_annotation/config.yaml`.


If everything is good and set, you can run the snakemake pipeline.
```bash
# Local 
snakemake -s mag_annotation/Snakefile --use-conda -j 1
```

Dependent on the number of samples, it may be relevant to run the Snake-flow on a High performance computing (HPC) server.
```bash
# HPC - this won't work unless you specify a legit group on your HPC in `config.yaml`
snakemake -s Snakefile --cluster qsub -j 32 --use-conda
```

### Workflow content
The workflow does the following. 
1. Splits the combined assembly into separate sample-specific ones 
2. Predicts proteins de novo using Prodigal meta
3. Searches proteins with VOGdb (PVOG) and miComplete bacterial hallmarks database 
4. Scans each contig using DeepVirFinder
5. Aggregates the above into Bin-wise annotation summaries used as input for Viral Decontamination Random Forest model

### Outputs
Running the provided script, the virome bins are written to a fasta file and bin-annotations are summarised in `vambbins_aggregated_annotation.txt`. Concatenate annotation files prior to. 
```bash
mkdir annotations
cat sample_annotation/*/*hmmMiComplete105.tbl > annotations/all.hmmMiComplete105.tbl
cat sample_annotation/*/*hmmVOG.tbl > annotations/all.hmmVOG.tbl
cat sample_annotation/*/*_dvf/*dvfpred.txt > annotations/all.DVF.predictions.txt
python mag_annotation/scripts/run_RF.py contigs.fna.gz vamb/clusters.tsv annotations resultdir
```
  
```bash
ls resultdir
resultdir/vambbins_aggregated_annotation.txt
resultdir/vambbins_RF_predictions.txt
resultsdir/vamb_bins #Concatenated putative viral bins - writes bins in chunks to files
```
 
We recommend VAMB bins to be evaluated with a dedicated Viral evaluation tool like [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) to identify HQ viruses.

```bash
checkv end_to_end resultsdir/vamb_bins/vamb_bins.1.fna checkv_vamb_bins  
```



### Further information 


- The RF model automates filtering of VAMB bins that are most likely bacterial and therefore provides a space of plausible viral entities for further validation. The RF-model has been trained on paired Metaviromes and Metagenomes to make precisde decisions based on simple parameteres as the ones below. Compared to a single contig viral prediction model, the RF approach is very accurate. The increased performance is likely explained by the RF model evaluating on bin-level where one sequence with a low viral score does NOT lead to a misprediction of the whole bin. Aggregated information (assuming the binning is really good!) from multiple-contigs improves prediction compared to single-contigs.

The RF model take few variables to make an accurate distinction.
| binsize (bp) | nhallm | distinct_VOGs_factor | cluster_DVF_score |
|--------------|--------|-------|-------------------|
| 2.000.000    | 100    | 0.2   | 0.3               |
| 60.000       | 3      | 1.3   | 0.7               |


- Bacterial MAGs and viral MAGs from the same metagenome can be efficiently associated using crispr-spacer approaches and sequence alignment (recommended cutoffs can be found in the article). From this, Host-viral abundance dynamics and bacterial pangenome modulation can be studied. Downstream viral proteome analysis should be based on the `viral regions` found in the `contamination.tsv` file produced by CheckV to prevent contaminating bacterial genes to influence the analysis. 
