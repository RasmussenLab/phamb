#!/bin/python
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='''
    - Parse the outfile from CheckM 
    - Writes out subset of Bins with a completeness >= 50 perc and contamination <= 10 perc
   ''')
parser.add_argument('-c', help='checkm file')
parser.add_argument('-o', help='checkm subset file')


def clean_checkm(args):
    checkmfileout = args.o
    best_MAGS = dict()

    ### Parse Checkm File
    with open(args.c, 'r') as infile, open(checkmfileout,'w') as outfile:
        outfile.write('bin\tmarker_lineage\t#genomes\t#markers\t#marker_sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tstrain_heterogenity\n')
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
            cluster = binid.split('_')[1]

            if float(completeness) >= 50 and float(contamination) <= 10:
                    if not cluster in best_MAGS:
                       best_MAGS[cluster] = (completeness,contamination,line)
                    else:
                       if float(completeness) >= best_MAGS[cluster][0] and float(contamination) <= best_MAGS[cluster][1]:
                          best_MAGS[cluster] = (completeness,contamination,line)

                    outfile.write('\t'.join(line)+'\n')

if __name__ == "__main__":
    args = parser.parse_args()
    clean_checkm(args)
    print('Done\n')
