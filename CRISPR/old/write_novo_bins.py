#!/bin/python
import argparse
import os
import sys
from Bio import SeqIO
import subprocess
import numpy as np
import gzip


parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-t',help='Directry with CRISPR cluster matrices.')
parser.add_argument('-i',help='CheckV results directory')
parser.add_argument('-o',help='Out Directory for .fna files')



def get_MAG_overview(args):
    """[Read CheckM file of MQNC MAGs]

    Args:
        args ([type]): [description]
    """
    MAGs = dict()
    checkm_file = '07_binannotation/bacteria/checkm/HMP2.all.MQNC.MAGS'
    with open(checkm_file,'r') as infile:
        header = infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            genomeid = line[0]
            motherbin = genomeid.split('_')[1]
            completeness = float(line[-3])
            qual = None
            if completeness >= 90:
                qual = 'HQ'
            else:
                qual = 'MQ'
            MAGs[genomeid] = dict()

    return MAGs

def read_cluster_viruses(args,bacterial_clusters):
    '''
    Loads Viral - Bacterial association matrix based on CRISPR

    '''
    viral_motherbins = set(['653','825','5971','3415','24810'])
    CRISPR_mates = {k:[] for k in bacterial_clusters}
    for bacterial_motherbin in bacterial_clusters:
        matrixfile = os.path.join(args.t,bacterial_motherbin+'.clustermatrix.txt')

        with open(matrixfile,'r') as infile:
            infile.readline()
            for line in infile:
                viral_cluster = line.strip().split('\t')[0]
                hits = line.strip().split('\t')[1:]
                hits = np.array(hits,dtype='int64')
                number_of_hits =  hits.sum()
                if number_of_hits > 1:
                    viral_motherbins.add(viral_cluster)
                    CRISPR_mates[bacterial_motherbin] += [viral_cluster]
    
    ### Parse CheckV standard file 
    quality_bins = {k:set() for k in viral_motherbins}
    checkv_quality = os.path.join(args.i,'quality_summary.tsv')
    with open(checkv_quality,'r') as infile:
        header = infile.readline().strip().split('\t')
        genome_copies_index = header.index('genome_copies')
        contamination_index = header.index('contamination')
        quality_index = header.index('checkv_quality')
        method_index = header.index('completeness_method')
        miuvig_quality = header.index('miuvig_quality')
        viral_genes_index = header.index('viral_genes')
        prophage_index = header.index('prophage')
        completeness_index = header.index('completeness')
        
        for line in infile:
            line = line.strip().split('\t')
            binid = line[0]
            motherbin = binid.split('_')[1]

            if not motherbin in viral_motherbins:
                continue
                
            quality = line[quality_index]
            completeness = line[completeness_index]
            genome_copies = float(line[genome_copies_index])
            contaminaton = line[contamination_index]
            

            ### Filter off Bins with little Viral evidence
            if quality == 'Not-determined' and line[viral_genes_index] == '0':
                    continue
            if completeness == 'NA':
                    continue
            if float(completeness) >= 50 and genome_copies < 1.25:
                    quality_bins[motherbin].add(binid)
    return CRISPR_mates,quality_bins

def write_cluster_fasta(args,quality_bins):
    '''
    Parse cleaned_contiggs.fna and write HQ-Bins to .fna file 
    Also write files with Binpaths
    
    '''

    ### Load Entries from cleaned_contigs.fna 
    bins_fasta = os.path.join(args.i,'cleaned_contigs.fna')

    bins_dict = {}
    for record in SeqIO.parse(open(bins_fasta, 'r'), 'fasta'):
        recordname = record.id
        motherbin = recordname.split('_')[1]
        bin = recordname
        if motherbin in quality_bins:
                if recordname in quality_bins[motherbin]:
                    if not motherbin in bins_dict:
                        bins_dict[motherbin] = {}
                        bins_dict[motherbin][bin] = record
                    else:
                        bins_dict[motherbin][bin] = record

    ### Write out Entries to seperate .fna files
    if not os.path.exists(args.o): 
        os.mkdir(args.o) 
    
    for motherbin in bins_dict:
        viral_bins_file = os.path.join('viral_cluster_files','vir'+motherbin+'.genomes.txt')

        with open(viral_bins_file,'w') as outfile:
            clusterdirectory = os.path.join(args.o,motherbin)
            if not os.path.exists(clusterdirectory):
                    os.mkdir(clusterdirectory)
            for bin in bins_dict[motherbin]:
                    binfile = os.path.join(clusterdirectory,bin+'.fna')
                    with open(binfile, 'w') as outhandle:
                        record = bins_dict[motherbin][bin]
                        SeqIO.write(record, outhandle, 'fasta')
                    outfile.write(binfile+'\n')



def prepare_bacterial_genome_lists(args,bacterial_clusters):
    """[summary]

    Args:
        args ([type]): [description]
        bacterial_clusters ([type]): [description]
    """

    MAGs = dict()
    checkm_file = 'bacteria/checkm/HMP2.all.MQNC.MAGS'
    with open(checkm_file,'r') as infile:
        header = infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            genomeid = line[0]
            motherbin = genomeid.split('_')[1]
            completeness = float(line[-3])
            contamination = float(line[-2])
            qual = None
            if completeness >= 50 and contamination <=20:
                MAGs[genomeid] = None

    if not os.path.exists('bacclusters'):
        os.system('mkdir -p bacclusters')

    bacterial_genomes = [f for f in os.listdir('bacteria') if '.fna' in f]
    bacterial_genomes_dict = dict()
    for genome in bacterial_genomes:

        motherbin = genome.replace('.fna','').split('_')[1]

        if not motherbin in bacterial_genomes_dict:
            bacterial_genomes_dict[motherbin] = [genome]
        else:
            bacterial_genomes_dict[motherbin] += [genome]

    for MAG in bacterial_clusters:
        fileout = os.path.join('bacclusters','bac'+MAG+'.genomes.txt')
        with open(fileout,'w') as outfile:
            for genome in bacterial_genomes_dict[MAG]:
                bin = genome.replace('.fna','')
                if bin in MAGs:
                    lineout = os.path.join('bacteria',genome)
                    outfile.write(lineout+'\n')


def run_cdhit(infile,fileout):
        try:
            executable = '/services/tools/cd-hit/4.8.1/cd-hit-est'
            command = [executable,
                    '-i',infile,
                    '-o',fileout,
                    '-d','0',
                    '-M','0',
                    '-aS','0.9']
            subprocess.check_call(command)
        except:
            message = 'ERROR: CD-hit finished abnormally.'
            print(message)
            sys.exit(1)


def run_fastani(bacterial_clusters,viral_cluster):
    
    executable='/services/tools/fastani/1.1/fastANI'
    bacteria_bins_file = os.path.join('bacclusters','bac'+bacterial_clusters+'.genomes.txt')
    if not os.path.exists(bacteria_bins_file):
        print('Bacterial genomes of cluster:',bacterial_clusters,'not there? Stopping this.')
        sys.exit(1)

    ### Make a temporary file with all Viral bins associated to the bacterial cluster
    viral_bins = []
    for viral_cluster in CRISPR_mates[bacterial_clusters]:
        viral_cluster_file = os.path.join('viral_cluster_files','vir'+viral_cluster+'.genomes.txt') 
        if os.path.exists( viral_cluster_file):
            with open(viral_cluster_file,'r') as infile:
                for line in infile:
                    line = line.strip()
                    viral_bins.append(line)

    viral_bins_file = os.path.join('viral_cluster_files',bacterial_clusters+'.allviruses.txt') 
    with open(viral_bins_file ,'w') as outfile:
        for viralbin in viral_bins:
            outfile.write(viralbin+'\n')
    
    print('Running fastani with:',bacterial_clusters)
    outfile = os.path.join('fastani','bac'+bacterial_clusters+'.viruses.fastani.txt')
    if not os.path.exists(outfile):
        if not os.path.exists('fastani'):
            os.system('mkdir -p fastani')
        try:
            command = [executable,
                    '-t','20',
                    '--fragLen','500',
                    '--minFrag','2',
                    '--ql',viral_bins_file,
                    '--rl',bacteria_bins_file,
                    '-o',outfile]
            subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        except:
            message = 'ERROR: fastaANI finished abnormally.'
            print(message)
            sys.exit(1)




def run_blastn(bacterial_clusters):
    """[Make FASTA files of all Genomes belonging to the given bacterial clcustser - run blastN towards all CRISPR associated virusses]

    Args:
        bacterial_clusters ([type]): [description]
    """

    
    bacteria_bins_file = os.path.join('bacclusters','bac'+bacterial_clusters+'.genomes.txt')
    
    ### Write Combined fasta for the Bacterial Cluster
    bacfiles = []
    with open(bacteria_bins_file,'r') as infile:
        for line in infile:
            bacfiles.append(line.strip())
    bacterial_genomes = os.path.join('bacclusters',bacterial_clusters+'.genomes.fna')
    bacterial_genomes_mapping = os.path.join('bacclusters',bacterial_clusters+'.contig_to_bin.txt')
    with open(bacterial_genomes, 'w') as outhandle, open(bacterial_genomes_mapping,'w') as outfile:
        for fastafile in bacfiles:
            binname = os.path.basename(fastafile).replace('.fna','')
            for record in SeqIO.parse(open(fastafile, 'r'), 'fasta'):
                header = record.id
                SeqIO.write(record, outhandle, 'fasta')
                outfile.write(header+'\t'+binname+'\n')
    
    ### Write fasta for Viral cluster
    viral_bins_file = os.path.join('viral_cluster_files',bacterial_clusters+'.allviruses.txt')
    viralfiles = []
    with open(viral_bins_file,'r') as infile:
        for line in infile:
            viralfiles.append(line.strip())
    viral_genomes = os.path.join('viral_cluster_files',bacterial_clusters+'.allviruses.fna')
    with open(viral_genomes, 'w') as outhandle:
        for fastafile in viralfiles:
            for record in SeqIO.parse(open(fastafile, 'r'), 'fasta'):
                SeqIO.write(record, outhandle, 'fasta')

    ### Index Bacterial genome fasta file
    outfile = os.path.join('blastn',bacterial_clusters+'.m6')

    if not os.path.exists(outfile):
        os.system('/services/tools/ncbi-blast/2.8.1+/bin/makeblastdb -in {} -dbtype nucl'.format(bacterial_genomes))

    executable='/services/tools/ncbi-blast/2.8.1+/bin/blastn'
    if not os.path.exists(outfile):
        if not os.path.exists('blastn'):
            os.system('mkdir -p blastn')
        try:
            command = [executable,
                    '-task','megablast',
                    '-evalue','0.001',
                    '-perc_identity','95',
                    '-db',bacterial_genomes,
                    '-out',outfile,
                    '-query',viral_genomes,
                    '-outfmt','6 std qlen slen',
                    '-num_threads','20']
            subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        except:
            message = 'ERROR: blastn finished abnormally.'
            print(message)
            sys.exit(1)
    
    ### Parse Blast results 
    if os.path.exists(outfile):
        parse_blastn_file(bacterial_clusters)

    ### Run CD-Hit To get Contig clusters
    #cdhit_fileout = os.path.join('bacclusters',bacterial_clusters+'.genomes.cdhit')
    #run_cdhit(bacterial_genomes,cdhit_fileout)



def run_blastn_prophages(bacterial_clusters,viral_prophages):
    """[Make FASTA files of all Genomes belonging to the given bacterial clcustser - run blastN towards all CRISPR associated virusses]

    Args:
        bacterial_clusters ([type]): [description]
    """

    
    bacteria_bins_file = os.path.join('bacclusters','bac'+bacterial_clusters+'.genomes.txt')
    
    ### Write Combined fasta for the Bacterial Cluster
    bacfiles = []
    with open(bacteria_bins_file,'r') as infile:
        for line in infile:
            bacfiles.append(line.strip())

    bacterial_genomes = os.path.join('bacclusters',bacterial_clusters+'.genomes.fna')
    bacterial_genomes_mapping = os.path.join('bacclusters',bacterial_clusters+'.contig_to_bin.txt')
    with open(bacterial_genomes, 'w') as outhandle, open(bacterial_genomes_mapping,'w') as outfile:
        for fastafile in bacfiles:
            binname = os.path.basename(fastafile).replace('.fna','')
            for record in SeqIO.parse(open(fastafile, 'r'), 'fasta'):
                header = record.id
                if len(record) > 10000: 
                    SeqIO.write(record, outhandle, 'fasta')
                    outfile.write(header+'\t'+binname+'\n')
    
    ### Write fasta for Viral cluster
    viral_bins_file = os.path.join('viral_cluster_files',bacterial_clusters+'.allviruses.txt')
    viralfiles = []
    with open(viral_bins_file,'r') as infile:
        for line in infile:
            viralfiles.append(line.strip())
    viral_genomes = os.path.join('viral_cluster_files',bacterial_clusters+'.allviruses.fna')
    with open(viral_genomes, 'w') as outhandle:
        for fastafile in viralfiles:
            for record in SeqIO.parse(open(fastafile, 'r'), 'fasta'):
                SeqIO.write(record, outhandle, 'fasta')

    ### Index Bacterial genome fasta file
    outfile = os.path.join('blastn',bacterial_clusters+'.m6')

    if not os.path.exists(outfile):
        os.system('/services/tools/ncbi-blast/2.8.1+/bin/makeblastdb -in {} -dbtype nucl'.format(bacterial_genomes))

    executable='/services/tools/ncbi-blast/2.8.1+/bin/blastn'
    if not os.path.exists(outfile):
        if not os.path.exists('blastn'):
            os.system('mkdir -p blastn')
        try:
            command = [executable,
                    '-task','megablast',
                    '-evalue','0.001',
                    '-perc_identity','95',
                    '-db',bacterial_genomes,
                    '-out',outfile,
                    '-query',viral_genomes,
                    '-outfmt','6 std qlen slen',
                    '-num_threads','20']
            subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        except:
            message = 'ERROR: blastn finished abnormally.'
            print(message)
            sys.exit(1)
    
    ### Parse Blast results 
    if os.path.exists(outfile):
        parse_blastn_file(bacterial_clusters)

    ### Run CD-Hit To get Contig clusters
    cdhit_fileout = os.path.join('bacclusters',bacterial_clusters+'.genomes.cdhit')
    run_cdhit(bacterial_genomes,cdhit_fileout)




def parse_blastn_file(bacterial_clusters):
    
    bacterial_genomes_mapping = os.path.join('bacclusters',bacterial_clusters+'.contig_to_bin.txt')
    contig_bin = dict()
    with open(bacterial_genomes_mapping,'r') as infile:
        for line in infile:
            contig, binid = line.strip().split('\t')
            if not contig in contig_bin:
                contig_bin[contig] = binid

    ### Parse Blast-results
    blastfile = os.path.join('blastn',bacterial_clusters+'.m6')
    Viral_genome_coverage = dict()
    with open(blastfile,'r') as infile:
        for line in infile:
            line = line.strip().split()
            viralbin = line[0]
            contig = line[1]
            alignment_length = int(line[3])
            
            viralbin_length = int(line[-2])
  
            if not viralbin in Viral_genome_coverage:
                Viral_genome_coverage[viralbin] = dict()
                Viral_genome_coverage[viralbin][contig] = np.zeros(viralbin_length)
            else:
                if not contig in Viral_genome_coverage[viralbin]:
                    Viral_genome_coverage[viralbin][contig] = np.zeros(viralbin_length)
    
            ref_from = int(line[8])
            ref_to = int(line[9])
            query_from = int(line[6])
            query_to = int(line[7])

            if query_to > query_from:
                Viral_genome_coverage[viralbin][contig][query_from:query_to] += 1
            else:
                Viral_genome_coverage[viralbin][contig][query_to:query_from] += 1
    
    blast_prophage_map = os.path.join('blastn',bacterial_clusters+'.prophage_in_MAG.txt') 
    with open(blast_prophage_map,'w') as outfile:
        for viralbin in Viral_genome_coverage:
            for contig in Viral_genome_coverage[viralbin]:
                bacterialbin = contig_bin[contig]
                bacterial_length = contig.split('_')[4]
                cov = Viral_genome_coverage[viralbin][contig]
                viralbin_length = len(cov)
                covered = 1-(len(cov[cov==0])/len(cov))
                if covered >= 0.5:
                    lineout = '\t'.join([viralbin, contig, bacterialbin, str(covered),str(viralbin_length), bacterial_length ])
                    outfile.write(lineout+'\n')
            






def parse_fastani_write_prophage(CRISPR_mates):
    """[summary]

    Args:
        CRISPR_mates ([type]): [description]

    Returns:
        [type]: [description]
    """

    if not os.path.exists('prophageannotations'):
        os.system('mkdir -p prophageannotations')


    viral_cluster_prophages = dict()
    for bacterial_clusters in CRISPR_mates:
        fileout = os.path.join('prophageannotations',bacterial_clusters+'.prophage.annotations.txt')
        with open(fileout,'w') as outfile:
            fastanifile = os.path.join('fastani','bac'+bacterial_clusters+'.viruses.fastani.txt')
            if os.path.exists(fastanifile):
                with open(fastanifile,'r') as infile:
                    for line in infile:
                        line = line.strip().split('\t')
                        ANI = float(line[2])
                        coverage = (int(line[3])/int(line[4]))*100

                        if ANI >= 95 and coverage >= 75:
                            binid = os.path.basename(line[0]).replace('.fna','')
                            viral_cluster = binid.split('_')[1]
                            viral_sample = binid.split('_')[0]
                            bacterialbin = os.path.basename(line[1]).replace('.fna','')

                            if not viral_cluster in viral_cluster_prophages:
                                viral_cluster_prophages[viral_cluster] = dict()
                                viral_cluster_prophages[viral_cluster][bacterialbin] = 1
                            else:
                                viral_cluster_prophages[viral_cluster][bacterialbin] = 1

                            sample = bacterialbin.split('_')[0]
                            lineout = '\t'.join([bacterial_clusters,bacterialbin,sample,viral_cluster,viral_sample,binid,str(ANI),str(coverage)])
                            outfile.write(lineout+'\n')
    return viral_cluster_prophages

def refine_clustermatrices(args,bacterial_clusters,viral_cluster_prophages):
    

    for bacterial_motherbin in bacterial_clusters:
        matrixfile = os.path.join(args.t,bacterial_motherbin+'.clustermatrix.txt')
        refined_matrixfile = os.path.join('prophageannotations',bacterial_motherbin+'.clustermatrix.prophage.txt')

        refined_matrix = []
        viral_clusters = []
        header = None
        with open(matrixfile,'r') as infile:
            header = infile.readline().strip()
            bacterialbins = header.split('\t')[1:]
            for line in infile:
                viral_cluster = line.strip().split('\t')[0]
                viral_clusters.append(viral_cluster)
                hits = line.strip().split('\t')[1:]
                hits = np.array(hits,dtype='float64')
                if viral_cluster in viral_cluster_prophages:
                    for i,bacterial_bin in enumerate(bacterialbins):
                        if bacterial_bin in viral_cluster_prophages[viral_cluster]:
                            hits[i] = hits[i] + 0.5
                hits = list(hits)
                refined_matrix.append(hits)
        
        refined_matrix = np.array(refined_matrix)
        refined_matrix = refined_matrix.astype('float64')

        viral_clusters = np.array(viral_clusters).reshape(-1,1)
        refined_matrix = np.column_stack((viral_clusters,refined_matrix))
        np.savetxt(refined_matrixfile, refined_matrix , delimiter='\t', header=header, newline='\n', fmt='%s')
                

    
def get_MAG_prophages(args):
    """[summary]

    Args:
        args ([type]): [description]
    """
    ### Parse MAG spacers
    MAGs = get_MAG_overview(args)


    ### Parse Prophage hits
    print('Parsing FastANI hits')
    prophage_file = os.path.join(args.p,'MGXVIR.fastani.txt.gz')
    with gzip.open(prophage_file,'rt') as infile:
        for line in infile:
            line = line.strip().split('\t')
            ANI = float(line[2])
            coverage = (int(line[3])/int(line[4]))*100
            if ANI >= 90 and coverage >= 75:
                phagebin = os.path.basename(line[0]).replace('.fna','')
                sample = phagebin.split('_')[0]
                phagemotherbin = phagebin.split('_')[1]
                MAG = os.path.basename(line[1]).replace('.fna','')
                bacterial_clusters = MAG.split('_')[1]

                if not phagemotherbin in MAGs[MAG]:
                    MAGs[MAG][phagemotherbin] = {'spacer':None,'prophage':True,'both':None}
                else:
                    if phagemotherbin in MAGs[MAG] and MAGs[MAG][phagemotherbin]['spacer'] is True:
                        MAGs[MAG][phagemotherbin] = {'spacer':True,'prophage':True,'both':True}


    



if __name__ == "__main__":
    args = parser.parse_args()

    #bacterial_clusters =['216','146']
    bacterial_clusters = [f.replace('.clustermatrix.txt','') for f in os.listdir(args.t)]
    #CRISPR_mates,quality_bins = read_cluster_viruses(args,bacterial_clusters)
    #write_cluster_fasta(args,quality_bins)
    
    ### Set up lists of bacterial genomes
    #prepare_bacterial_genome_lists(args,bacterial_clusters)
     

    ### Execute Fastani 
    #for motherbin in CRISPR_mates:
    #    print('FastANI for ',motherbin)
        #run_fastani(motherbin,CRISPR_mates)

    priorited_bacteria = ['146','216','177']
    for bac in priorited_bacteria:
        run_blastn(bac)

    ### Write Prophage Annotations for each Bacterial Motherbin 
    #viral_cluster_prophages = parse_fastani_write_prophage(CRISPR_mates)
    #refine_clustermatrices(args,bacterial_clusters,viral_cluster_prophages)

