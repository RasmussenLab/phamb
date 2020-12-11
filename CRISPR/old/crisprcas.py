#!/bin/python

import multiprocessing as mp
import subprocess
import numpy as np
import os
import sys
import json
import argparse
import time as _time


parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
   ''')
parser.add_argument('-b', help='Bacterial annotation directory')
parser.add_argument('-i', help='Directory with spacer-blast files')
parser.add_argument('-c', help='Checkm file')
parser.add_argument('-t', help='# Threads')

def run_crisprcas(args,binid):
    '''
    Parallelise crisprcas runs on MQ and HQ bins
    '''

    print('Parsing:{}'.format(binid))
    bindir = args.b
    try:
        executable = '/home/projects/cpr_10006/people/joacj/miniconda3/envs/crisprcasfinder/bin/perl'
        directory_out = os.path.join(bindir,'crisprcasfinder/ccf_'+binid,'ccf')
        os.system('mkdir -p {}'.format( os.path.join(bindir,'crisprcasfinder/ccf_'+binid ))  )
        if os.path.exists(os.path.join(directory_out,'result.json')): 
            ### Write out spacers to file: ccf/spacers.fna
            write_out_spacers(directory_out,binid)
        else:
            _time.sleep(1)
            command = [executable,
            os.path.join(bindir,'CRISPRCasFinder/CRISPRCasFinder-release-4.2.20/CRISPRCasFinder.pl'),
            '-out',directory_out,
            '-in', os.path.join(bindir, binid +'.fna'),
            '-quiet',
            '-fast',
            '-soFile','/home/people/joacj/bin/selv392v2.so',
            '-minSP','22']
            run_command = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            write_out_spacers(directory_out,binid)
    except:
        message = 'ERROR: Something went wrong'
        print(run_command)
        print(message)
        sys.exit(1)

def run_CRISPRDetect(args,binid):
    '''
    Parallelise crisprcas runs on MQ and HQ bins
    '''
    print('Parsing:{}'.format(binid))
    bindir = args.b
    try:
        executable = '/home/projects/cpr_10006/people/joacj/miniconda3/envs/crisprcasfinder/bin/perl'
        bindir2 = os.path.join(bindir,'crisprdetect/crisprdetect_'+binid)
        os.system('mkdir -p {}'.format(bindir2))
        os.chdir(bindir2)
        if os.path.exists(binid+'_CRISPRDetect'):
            print('CRISPRDetect file exists - writing spacers')
            ### Write out spacers to file: ccf/spacers.fna
            pass
        else:
            command = [executable,
            os.path.join(bindir,'CRISPRDetect.pl'),
            '-o',binid+'_CRISPRDetect',
            '-f',os.path.join(bindir,binid+'.fna'),
            '-check_direction','0',
            '-array_quality_score_cutoff','3',
            '-T','1']
            subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    except:
        message = 'ERROR: Something went wrong'
        print(message)
        sys.exit(1)




def write_out_spacers(directory_out,binid):
    with open(os.path.join(directory_out,'result.json')) as json_file:
        data = json.load(json_file)
    crispr_contigs = []
    for i,entry in enumerate(data['Sequences']):
        if len(entry['Crisprs']) >0:
            crispr_contigs.append(i) 
    with open(os.path.join(directory_out,'new.spacers.fna'),'w') as out:
        for i in crispr_contigs:
            crisp = data['Sequences'][i]
            contigid = crisp['Id']
            for cr in crisp['Crisprs']:
                evidence_level = cr['Evidence_Level']
                for j, reg in enumerate(cr['Regions']):
                    if reg['Type'] == 'Spacer':
                        spacerid = binid + ':' + str(j) + '_' + contigid + '_evidence_'+str(evidence_level)
                        out.write('>{}\n{}\n'.format(spacerid,reg['Sequence']))
    



def blastn():
    if os.path.exists('07_binannotation/bacteria/crisprcasfinder/new.all.spacers.fna'):
        try:
            executable = '/home/projects/cpr_10006/people/joacj/miniconda3/envs/crisprcasfinder/bin/blastn'
            db = '04_annotation/annotation_summaries/VAMBV3.Viral_RF_predictions.bins.fna'
            command = [executable,
                    '-task','blastn-short',
                    '-evalue','0.001',
                    '-perc_identity', '90',
                    '-num_threads','24',
                    '-db',db,
                    '-query','07_binannotation/bacteria/crisprcasfinder/new.all.spacers.fna',
                    '-out','07_binannotation/bacteria/crisprcasfinder/new.spacers.viralfragment.m6',
                    '-outfmt', '6 std qlen slen']
            subprocess.check_call(command)
        except:
            message = 'ERROR: BlastN finished abnormally.'
            print(message)
            sys.exit(1)

def cdhit(spacerfile,fileout):
        try:
            executable = '/services/tools/cd-hit/4.8.1/cd-hit-est'
            command = [executable,
                    '-i',spacerfile,
                    '-o',fileout,
                    '-d','0',
                    '-aS','0.9']
            subprocess.check_call(command)
        except:
            message = 'ERROR: CD-hit finished abnormally.'
            print(message)
            sys.exit(1)



def write_MAG_spacers(args,magspacerdict,spacers_in_MAG):
    """[summary]

    Args:
        args ([type]): [description]
        magspacerdict ([type]): [description]
    """

    spacer_dir = os.path.join(args.i,'spacer_by_MAG')
    if not os.path.exists(spacer_dir):
        os.system('mkdir -p {}'.format(spacer_dir))
    
    ### write out table with Metadata to each Spacer.
    metadata = {}
    metadata_file = 'metadata'
    if os.path.exists(metadata_file):
        with open(metadata_file,'r') as infile:
            header = infile.readline()
            for line in infile:
                project, Sample, Subject, week_num, diagnosis, datatype = line.strip().split('\t')
                metadata[Sample] = [Subject,week_num,diagnosis]
    
    for MAG in spacers_in_MAG:
        mag_spacer_dir =os.path.join(args.i,'spacer_by_MAG',MAG)
        if not os.path.exists(mag_spacer_dir):
            os.system('mkdir -p {}'.format(mag_spacer_dir))
        
        fileout = os.path.join(mag_spacer_dir,MAG+'.spacers.fna')
        with open(fileout,'w') as outfile:
            for spacerid in spacers_in_MAG[MAG]:
                seq = magspacerdict[spacerid]
                outfile.write('>'+spacerid +'\n' + seq + '\n')
        clstr_file = os.path.join(mag_spacer_dir,MAG+'.cdhit')
        if not os.path.exists(clstr_file):
            cdhit(fileout,clstr_file)
        
        ### Parse cluster annotation
        spacerclstr = dict()
        i = 0
        clstr = None
        clstr_file = os.path.join(mag_spacer_dir,MAG+'.cdhit.clstr')
        with open(clstr_file,'r') as infile:
            for line in infile:
                if line[0] == '>':
                    i += 1 
                    clstr = str(i)
                else:
                    line = line.strip().split('>')[1]
                    spacerid = line.split('...')[0]
                    spacerclstr[spacerid] = clstr
        fileout = os.path.join(mag_spacer_dir,MAG+'.spacers.annotated.fna')      
        with open(fileout,'w') as outfile:
            for spacerid in spacers_in_MAG[MAG]:
                clstrnumber = spacerclstr[spacerid]
                seq = magspacerdict[spacerid]
                Sample = spacerid.split(':')[0].split('_')[0]
                Subject,week_num,diagnosis = metadata[Sample]
                new_spacerid = Sample +'_'+ Subject + '_' + 'week:' +clstrnumber + '_' + 'spacergrp:' + clstrnumber
                outfile.write('>'+new_spacerid +'\n' + seq + '\n')


def write_out_arrays(args,binid):

    Sample = binid.split('_')[0]

    ### Read metadata
    metadata = {}
    metadata_file = 'metadata'
    if os.path.exists(metadata_file):
        with open(metadata_file,'r') as infile:
            header = infile.readline()
            for line in infile:
                project, Sample, Subject, week_num, diagnosis, datatype = line.strip().split('\t')
                metadata[Sample] = [Subject,week_num,diagnosis]


    Subject,week_num,diagnosis = metadata[Sample]

    crisprcasfile = os.path.join(args.i,'ccf'+'_'+binid,'ccf','result.json')
    casettefile = os.path.join(args.i,'ccf'+'_'+binid,'ccf','casette.fna')


    if os.path.exists(crisprcasfile):
        with open(crisprcasfile) as json_file:
            data = json.load(json_file)
        crispr_contigs = []
        for i,entry in enumerate(data['Sequences']):
            if len(entry['Crisprs']) >0:
                crispr_contigs.append(i)
        
        with open(casettefile,'w') as out:
            for i in crispr_contigs:
                crisp = data['Sequences'][i]
                contigid = crisp['Id']

                for cr in crisp['Crisprs']:
                    evidence_level = cr['Evidence_Level']
                    orientation = cr['Potential_Orientation']
                    arraystart = cr['Start']
                    arrayend = cr['End']
                    casette = ''
                    casette_header = Sample + '_' + Subject + '_' + week_num + '_' + contigid 
                    for j, reg in enumerate(cr['Regions']):
                        crisptype = reg['Type']
                        seq = reg['Sequence']
                        if crisptype == 'DR':
                            casette += '['+seq+']'
                        else:
                            casette += seq
                    out.write('>{}\n{}\n'.format(casette_header,casette))


def read_spacer_dict(args):
    """[Reads a .fna file with spacers and maps the spacer-sequence to the entry naame ]

    Args:
        args ([type]): [description]
    """
    spacerfile = os.path.join(args.i, 'new.all.spacers.fna')
    spacers_in_MAG = dict()
    magspacerdict = dict()
    spacerid = None
    with open(spacerfile,'r') as infile:
        for line in infile:
            if line[0] == '>':
                spacerid = line[1:].strip()
                bin = spacerid.split(':')[0]
                motherbin = bin.split('_')[1]
                if not motherbin in spacers_in_MAG:
                    spacers_in_MAG[motherbin] = [spacerid]
                else:
                    spacers_in_MAG[motherbin] += [spacerid]
                
            else:
                seq = line.strip()
                magspacerdict[spacerid] = seq
    write_MAG_spacers(args,magspacerdict,spacers_in_MAG)

    return magspacerdict



def write_MAG_matrix(args,selected_MAG,phage_MAG,MAGs):

    ### What's the dimensions?
    # Number of differeent Viral Motherbins x Number of Samples the MAG 
    mags = []
    genomes = list(MAGs[selected_MAG].keys())
    for binid in genomes:
        if binid == 'viral_clusters':
            continue
        hit = list(MAGs[selected_MAG][binid])
        if hit[0] == 1 and hit[1] == 'HQ':
            mags.append(binid)
    print(len(mags),' Bins in cluster targetting a Viral Cluster')

    if len(mags) > 1:
        
        ### Organise all MAGs associiated with the Cluster
        viral_cluster_mags = dict()
        for phagemotherbin in MAGs[selected_MAG]['viral_clusters']:
            crispr_mags = set()
            ### Get set of MAGs 
            for phagebin in phage_MAG[phagemotherbin]:
                for MAG in phage_MAG[phagemotherbin][phagebin]:
                    crispr_mags.add(MAG)
            viral_cluster_mags[phagemotherbin] = crispr_mags
        
        ### Now return binary hit-list 
        hitmatrix = []
        clusters = []
        for phagemotherbin in MAGs[selected_MAG]['viral_clusters']:
            crispr_mags = viral_cluster_mags[phagemotherbin]
            mag_hit_list = np.zeros(len(mags))
            for i,mag in enumerate(mags):
                if mag in crispr_mags:
                    mag_hit_list[i] = 1
            hitmatrix.append(mag_hit_list)
            clusters.append(phagemotherbin)
        
        header = ['mothercluster'] + mags
        header = '\t'.join(header)
        hitmatrix = np.array(hitmatrix)
        hitmatrix = hitmatrix.astype('int64')

        clusters = np.array(clusters).reshape(-1,1)
        hitmatrix = np.column_stack((clusters,hitmatrix))

        ### Write out matrix.
        os.system('mkdir -p 07_binannotation/bacteria/crisprcasfinder/clustermatices')
        tableout =  os.path.join('07_binannotation/bacteria/crisprcasfinder/clustermatices',selected_MAG + '.clustermatrix.' + 'txt')
        np.savetxt(tableout, hitmatrix , delimiter='\t', header=header, newline='\n', fmt='%s')



if __name__ == "__main__":
    args = parser.parse_args()


    ### Parse MQ and NC Bin's to process from CHECKM file
    bindir = args.b
    binids = []
    with open(args.c,'r') as infile:
        infile.readline()
        for line in infile:
            binid = line.strip().split('\t')[0]
            binids.append(binid)
    
    ### For each MAG, run CRISPR-spacer and array detection in Parallel
    ### Using Multiple threads! 
    processresults = list()
    if not os.path.exists(os.path.join(args.i,'all.spacers.fna')):
        print('Running CRISPRCASFINDER')
        n_cores = args.t

        with mp.Pool(processes=n_cores) as pool:
            for binid in binids:
                arguments = (args,binid)
                processresults.append(pool.apply(run_crisprcas, arguments))
            
            all_done, any_fail = False, False
            while not (all_done or any_fail):
                _time.sleep(5)
                all_done = all(process.ready() and process.successful() for process in processresults)
                any_fail = any(process.ready() and not process.successful() for process in processresults)
                if all_done:
                    pool.close()
            
            pool.join()

        ### Concatenate spacer files
        os.system('cat 07_binannotation/bacteria/crisprcasfinder/*/ccf/new.spacers.fna > 07_binannotation/bacteria/crisprcasfinder/all.spacers.fna')

    ### Write out arrays
    for binid in binids:
        casettefile = os.path.join(args.i,'ccf'+'_'+binid,'ccf','casette.fna')
        if not os.path.exists(casettefile):
            write_out_arrays(args,binid)

    ### Run Blastn
    if not os.path.exists(os.path.join(bindir,'crisprcasfinder/new.spacers.viralfragment.m6')):
        print('Running Blast')
        blastn()

