#!/bin/python

import os 
import sys
import subprocess
from Bio import SeqIO


###  Predict proteins of NC genomes
def write_out_proteomes_prodigal(args,ncbins,prodigalfile_AA):

    ### Parse checkV generated fasta with Host-contamination removed 
    checkv_directory = args.v

    ### Predict proteins of Created NC viral Genomes
    prodigalfile_GFF = os.path.join(checkv_directory,'nc_genomes_proteins.gff')
    ncgenomes_outfile = os.path.join(checkv_directory,'nc_genomes.fna')
    if os.path.exists(prodigalfile_AA):
        print('Proteins of NC genomes already predicted, delete first to recreate ', prodigalfile_AA)
    else:
        try:
            prodigalexecutable = '/services/tools/prodigal/2.6.3/prodigal'
            command = [prodigalexecutable,
                    '-i', ncgenomes_outfile,
                    '-a', prodigalfile_AA,
                    '-o', prodigalfile_GFF,
                    '-p', 'meta',
                    '-g', '11',
                    '-q',
                    '-f', 'gff']
            subprocess.check_call(command)
        except:
            message = 'ERROR: Prodigal finished abnormally.'
            print(message)
            sys.exit(1)


def write_out_proteomes(args,ncbins,prodigalfile_AA):

    ### Parse checkV generated fasta with Host-contamination removed 
    checkv_directory = args.v

    written_proteins = set()
    all_proteins = os.path.join(checkv_directory,'tmp/proteins.faa')
    outhandle = open(prodigalfile_AA , 'w')
    for record in SeqIO.parse(open(all_proteins, 'r'), 'fasta'):
        recordname = record.id
        binid = '_'.join(recordname.split('_')[:2])
        if binid in ncbins: 
            if recordname in written_proteins:
                continue
            written_proteins.add(recordname)
            record.description = ''
            SeqIO.write(record, outhandle, 'fasta')

###
def crass_blast(databasedirectory,executable, proteins, outfile):

    try:
        db = os.path.join(databasedirectory,'TaxonomyDB/RefSeq_2019-09-20/crassphage3/UGP_092.faa')
        command = [executable,
                '-task','blastp',
                '-evalue','0.0001',
                '-num_threads','19',
                '-db',db,
                '-query',proteins,
                '-out',outfile,
                '-outfmt', '6 std qlen slen']
        subprocess.check_call(command)
    except:
        message = 'ERROR: BlastP finished abnormally.'
        print(message)
        sys.exit(1)
    print('BlastP against Crassphage Proteins Complete.')



def PLSDB_blast(databasedirectory,executable,genomes,outfile):
    
    try:
        db = '/home/projects/cpr_10006/projects/phamb/tools/plsdb/data/master/2020_01_14.fna'
        command = [executable,
                '-task','megablast',
                '-evalue','0.0001',
                '-num_threads','19',
                '-perc_identity','85',
                '-db',db,
                '-query',genomes,
                '-out',outfile,
                '-outfmt', '6 std qlen slen']
        subprocess.check_call(command)
    except:
        message = 'ERROR: BlastN finished abnormally.'
        print(message)
        sys.exit(1)



### BlastP against Eukaryotic Viruses 

def herpes_blast(databasedirectory,executable, genomes, outfile):
    try:
        db = os.path.join(databasedirectory,'herpes_db/proteomes/nr.vir.euk')
        command = [executable,
                'blastx',
                '-d',db,
                '--sensitive',
                '--index-chunks','1',
                '--threads','19',
                '-q',genomes,
                '-o',outfile,
                '--outfmt', '6']
        
        subprocess.check_call(command)
    except:
        message = 'ERROR: Diamond BlastP finished abnormally.'
        print(message)
        sys.exit(1)
    print(' Diamond BlastP against Eukaryotic Viruses Complete.')


### Search against RVDB HMM profiles
def RVDB_search(databasedirectory,executable, proteins, outfile):
    try:
        db = os.path.join(databasedirectory,'rvdb/U-RVDBv18.0-prot.hmm')
        command = [executable,
                '-E','0.00001',
                '--cpu','19',
                '--tblout',outfile,
                db,
                proteins]
        subprocess.check_call(command,stdout=subprocess.DEVNULL)
    except:
        message = 'ERROR: Hmmsearch finished abnormally.'
        print(command)
        print(message)
        sys.exit(1)
    
    print('HMMsearch against RVDB database finished.')

### Search against VOG HMM profiles
# Based on this search - it will be much easier to write out Hallmark proteins for Phylogenetic Trees

def VOG_search(databasedirectory,executable, proteins, outfile):
    try:
        db = os.path.join(databasedirectory,'VOG/vog_hmms/AllVOG.hmm')
        command = [executable,
                '-E','0.00001',
                '--cpu','19',
                '--tblout',outfile,
                db,
                proteins]
        subprocess.check_call(command,stdout=subprocess.DEVNULL)
    except:
        message = 'ERROR: Hmmsearch finished abnormally.'
        print(message)
        sys.exit(1)
    
    print('HMMsearch against VOG database finished.')


### Usearch against UNIprot Viral Proteins
def uniprot_search(databasedirectory,executable, proteins, outfile):
    try:
        db = os.path.join(databasedirectory,'Demovir/uniprot_trembl.viral.udb')
        command = [executable,
                    '-db',db,
                    '-ublast',proteins,
                    '-evalue', '1e-5',
                    '-blast6out',outfile,
                    '-threads','19']
        subprocess.check_call(command)
    except:
        message = 'ERROR: Usearch finished abnormally.'
        print(message)
        sys.exit(1)
    print('Usearch against Viral Trembl Complete.')