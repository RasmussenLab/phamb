#!/bin/bash

### Usage 
# myqsub -d `pwd` -e log/cdhitest.e -o log/cdhitest.e 
#PBS -l nodes=1:ppn=20:thinnode
#PBS -l walltime=1:00:00:00
#PBS -l mem=80gb

module load perl/5.24.0
module load ncbi-blast/2.8.1+
module load cd-hit/4.8.1

cd-hit-est -i 07_binannotation/checkv/VAMB_bins/cleaned_contigs.fna -o 07_binannotation/checkv/VAMB_bins/cleaned_contigs.NR.fna -sc -aS 0.9 -d 0 -T 20 -M 0
