#!/bin/python 
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='''
    Reformatting fastani all vs all search of VAMB bins
''')

parser.add_argument('-f', help='fastani file')
parser.add_argument('-o', help='Formatted fastani file out')

def parse_format_fastani(fastanifile, outfile):
    with open(fastanifile,'r') as infile, open(outfile,'w') as out:
        outheader = ['query','reference','qcluster','rcluster','ani','bifrag','totalfrag','fastani_coverage']
        out.write('\t'.join(outheader)+'\n')
        for line in infile:
            query, reference, ani, bifragment, totalfragments = line.strip().split()
            query_coverage = round(int(bifragment)/int(totalfragments),2)            
            if float(ani) >=85:
                query = os.path.basename(query).replace('.fna','')                
                reference = os.path.basename(reference).replace('.fna','')    
                query_cluster = query.split('_')[1]
                reference_cluster = reference.split('_')[1]
                lineout = [query, reference, query_cluster,reference_cluster , ani, bifragment, totalfragments, str(query_coverage)] 
                out.write('\t'.join(lineout) + '\n')

if __name__ == "__main__":
    args = parser.parse_args()
    parse_format_fastani(fastanifile= args.f, outfile= args.o)