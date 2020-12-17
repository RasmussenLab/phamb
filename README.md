# phamb
A Phage discovery approach through binning of Metagenomic derived contigs 

## Prerequisites - Snakemake 

In order to run most of it, you need conda installed in order to install snakemake.

```
conda install -n snakemake snakemake pygraphviz python=3.8

```


## reads to bins - Snakemake pipeline 
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


## MAG annotation 
NB: Requires VAMB bins and assemblies etc. from above 

- Annotation of proteins from assembled contigs, including:
    - PVOG 
    - DeepVirFinder score for each contig
    - Bacterial Hallmark annotation
- Bin-wise annotation summaries used as input for Viral Decontamination

i.e.
Using relatively few variables automates filtering of contigs belonging to metagenomic bins that are NOT likely bacterial and therefore plausible viral/others.


| binsize (bp) | nhallm | nVOGs | cluster_DVF_score |
|--------------|--------|-------|-------------------|
| 2.000.000    | 100    | 0.2   | 0.3               |
| 60.000       | 3      | 1.3   | 0.7               |

## CRISPR
- Using CRISPR-cas-typer (CCtyper) to extract CRISPR-arrays from MAGs generated using VAMB
- Alignment of putative Viruses to MAGs
- Summarise Viral-MAG connections

## CHECKV (Work in progress)
Code for parsing checkV output files and annotating MQ/NC viral bins with taxonomical and functional annotation

## MAG VMAG Abundance 
Code for making abundance matrices of MAGs and Viruses (validated with CheckV)


