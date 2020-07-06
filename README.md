# phamb
Downstream processing of VAMB binning for Viral Elucidation


## reads to bins 
A series of pipeline steps written in bash that runs the following on a metagenomic dataset of choice
- QC
- Assembly
- Mapping
- Annotation 
Pipeline options ara available in config.

All steps necessary for running metagenomic binning and the annotation provides some foundation for running decontamination of Viral contigs. 

## Viral decon(tamination)
Using relatively few variables automates filtering of contigs belonging to metagenomic bins that are NOT likely bacterial and therefore plausible viral/others.

i.e.

| binsize (bp) | nhallm | nVOGs | cluster_DVF_score |
|--------------|--------|-------|-------------------|
| 2.000.000    | 100    | 0.2   | 0.3               |
| 60.000       | 3      | 1.3   | 0.7               |

## CHECKV 
Code for parsing checkv output files and annotating MQ/NC viral bins with taxonomical and functional annotation



