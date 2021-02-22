# phamb
A Phage discovery approach through binning of Metagenomic derived contigs. Most snakemake workflows comes with conda-environments, thus dependencies and programmes are automatcally installed. 

## Prerequisites - Snakemake 

In order to run most of it, you need conda installed in order to install snakemake.

```
conda install -n snakemake snakemake pygraphviz python=3.8

```


## MAG annotation for isolating Viral bins  

### Requirements
VAMB bins and concatenated assemblies. 

```
contigs.fna.gz #Concatenated assembly 
vamb/clusters.tsv   #Clustered contigs based on the above contigs.fna.gz file 
```

Furthermore. 
* VOGdb (https://vogdb.csb.univie.ac.at/download) - untar `vog.hmm.tar.gz` to get `AllVOG.hmm`    [File needs to be specified in `config.yaml`]
* Micomplete Bacterial HMMs (https://bitbucket.org/evolegiolab/micomplete/src/master/micomplete/share/Bact105.hmm)   [File needs to be specified in `config.yaml`]
* Clone DeepVirFinder (git clone https://github.com/jessieren/DeepVirFinder) 


### How to Run 

```
mkdir -p projectdir 
cd projectdir 
git clone https://github.com/RasmussenLab/phamb.git
cp -r phamb/workflows/mag_annotation .
python mag_annotation/scripts/split_contigs.py -c contigs.fna.gz 

```

Now the `contigs.fna.gz` is splitted into individual assemblies i.e. `assembly/{sample}/{sample}.fna`
In addition, a `sample_table.txt` file is created with a line for each sample.
Check that `sample_table.txt` contains sample identifiers corresponding to the ones you expect. 
The number of lines should correspond to the number of samples used to make the concatenated assembly

Now, Specify paths for databases, vamb directory, location of assembly  and computational resouces in `mag_annotation/config.yaml`  


If everything good and set, you can run the snakemake pipeline.
```
# Local 
snakemake -s mag_annotation/Snakefile --use-conda 
```

```
# HPC - this won't work unless you specify a legit group on your HPC in `config.yaml`
snakemake -s Snakefile --cluster qsub -j 32 --use-conda 

```


### Workflow content
The workflow does the following. 
1. Splits the combined assembly into separate sample-specific ones 
2. Predicts proteins de novo using Prodigal meta
3. Searches proteins with VOGdb (PVOG) and miComplete bacterial hallmarks db 
4. Scans each contig using DeepVirFinder  
5. Aggregates the above into Bin-wise annotation summaries used as input for Viral Decontamination Random Forrest model

By the end of the day, what goes on is pretty straight forward. Using relatively few variables, the RF automates filtering of VAMB bins that are most likely bacterial and there provides a space of plausible viral/plasmid entities for further validation. The RF-model has been trained on paired Metaviromes and Metagenomes to make precisde decisions based on simple parameteres as the ones below. 

| binsize (bp) | nhallm | nVOGs | cluster_DVF_score |
|--------------|--------|-------|-------------------|
| 2.000.000    | 100    | 0.2   | 0.3               |
| 60.000       | 3      | 1.3   | 0.7               |

Why does this work so well? Aggregated information (assuming the binning is really good!) from multiple-contigs simplifies prediction compared to single-contigs.


## reads to bins - If you are starting from Scratch with metagenomes or metaviromes 
A series of pipeline steps that runs the following on a metagenomic dataset of choice
- QC
- Assembly
- Mapping
- Binning using VAMB 

Pipeline options ara available in config.
Most necessary requirements are packed into Conda-environments, some are not currently.

```
snakemake -s crispr/Snakefile -j --use-conda --use-envmodules
```

## CRISPR
- Using CRISPR-cas-typer (CCtyper) to extract CRISPR-arrays from MAGs generated using VAMB
- Alignment of putative Viruses to MAGs
- Summarise Viral-MAG connections

## CHECKV (Work in progress)
Code for parsing checkV output files and annotating MQ/NC viral bins with taxonomical and functional annotation

## MAG VMAG Abundance 
Code for making abundance matrices of MAGs and Viruses (validated with CheckV)


