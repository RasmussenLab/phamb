#!/usr/bin/env python
import os
import sys 
'''
Parses a fasta, splits it into individual .fna files with identifier as name.
'''
infile = sys.argv[1]
outdir = sys.argv[2]

f=infile
if os.path.exists(f):
        opened=False
        with open(f,'r') as infile:
                for line in infile:
                    if line[0] == ">":
                        if opened:
                            of.close()
                        opened = True
                        header = line[1:].strip()
                        of=open("{}/{}.fna".format(outdir,header), "w")
                    of.write(line)
        of.close()