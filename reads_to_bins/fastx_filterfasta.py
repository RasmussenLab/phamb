#!/usr/bin/env python

from __future__ import division

import argparse
import os
import logging
import re
from Bio import SeqIO
import sys

# create the parser
parser = argparse.ArgumentParser(description='''
   Filter fa file based on sequence lengths
   ''')

# add the arguments
parser.add_argument('--i', help='file to be filtered')
parser.add_argument('--o', help='outfile')
parser.add_argument('--min', help='minimum length', type=int, default=100)
parser.add_argument('--id', help='sample id',type=str,default=None)

# parse the command line
args = parser.parse_args()
#args = parser.parse_args('--i graph_prefix_k29.contig --min 100'.split())

# set paths
home = os.getcwd()

# get filenames
file = args.i
path = os.path.split(file)[0]
filename = os.path.split(file)[1]
fileprefix = os.path.splitext(filename)[0]
pattern = re.compile(r'(.+?)(\.[^.]*$|$)')
match = pattern.search(fileprefix)
if match:
   basename = match.group(1)
else:
   basename = fileprefix

# parse
outfile = args.o
outhandle = open(outfile , 'w')

sample_name = args.id
pass_count = 0
total_count = 0
for record in SeqIO.parse(open(args.i, 'r'), 'fasta'):
   total_count += 1
   if len(record) >= args.min:
      new_record = sample_name + '_' + record.id + '_contig_' + str(pass_count) 
      record.id = new_record 
      record.description = ''
      SeqIO.write(record, outhandle, 'fasta')
      pass_count += 1

sys.stderr.write('%s: Written %i of %i (%f)\n' % (outfile, pass_count, total_count, pass_count/total_count*100))

