#!/bin/python
import argparse
import os
import sys
import pathlib
import taxopy
from Bio import SeqIO
import subprocess
import re
import numpy as np

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-c',help='VAMB clusterfile/clusters.tsv')
parser.add_argument('-a',help='Annotation file Direcotry')
parser.add_argument('-v', help='checkv directory')



def organise_functional_annotations(args,ncbins):

    def load_bin_proteins(args):
        
        proteinfile = os.path.join(args.v,'nc_genomes_proteins.faa')

        bin_proteomes = {}
        with open(proteinfile,'r') as infile:
            for line in infile:
                if line[0] =='>':
                    proteinid = line.split('#')[0].strip().replace('>','')
                    binid = '_'.join(proteinid.split('_')[:2])
                    

    
    def load_KO_annotation():
        '''Load description for each KO'''
        kodbfile = os.path.join('/home/projects/cpr_10006/projects/phamb/databases/kofam','ko_list')

        kodict = {}
        with open(kodbfile,'r') as infile:
            header = infile.readline().strip().split('\t')
            knum_index = header.index('knum')
            definition_index = header.index('definition')
            for line in infile:
                line = line.strip().split('\t')
                knum = line[knum_index]
                definition = line[definition_index]
                kodict[knum] = definition
        return kodict
    
    def load_GO_mapping():
        '''Load GO-slim information for each GO. Plus, get the GO-slim annotation from interpro2go flat-file'''
        goslim = {}
        goterm, name, namespace = None, None, None
        dbdir = '/home/projects/cpr_10006/projects/phamb/databases/interpro'
        goslim_files = ['goslim_metagenomics.obo','goslim_agr.obo']
        for gofile in goslim_files:
            with open(os.path.join(dbdir,gofile),'r') as infile:
                for line in infile:
                    if line[:3] == 'id:':
                        goterm = line.strip().split(' ')[1]
                        goslim[goterm] = {'name':None,'namespace':None}
                    elif line[:5] == 'name:':
                        name = line.strip().split(':')[1]
                        name = name.strip()
                        goslim[goterm]['name'] = name
                    elif line[:9] == 'namespace':
                        namespace = line.strip().split(':')[1]
                        namespace = namespace.strip()
                        goslim[goterm]['namespace'] = namespace
                    if line[:6] == 'alt_id':
                        altgoterm = line.strip().split(' ')[1]
                        goslim[altgoterm] = {'name':None,'namespace':None}
                        goslim[altgoterm]['name'] = name
                        goslim[altgoterm]['namespace'] = namespace
        goslim = {k:v for k,v in goslim.items() if 'GO' in k}
        
        ### Parse interpro2go annotation
        interpro2go = {}
        with open(os.path.join(dbdir,'interpro2go')) as infile:
            for line in infile:
                if line[0] == '!':
                    continue       
                line = line.strip()
                iprid= line.split('>')[0].split(' ')[0]
                iprid = iprid.replace('InterPro:','')
                dircount = line.count('>')
                if dircount == 1:
                    high_level_function, GOid = line.split('>')[1].split(';')
                    high_level_function = high_level_function.replace('GO:','').strip()
                    interpro2go[iprid] = {'GO':GOid, 'GO_low':high_level_function}
                else:
                    print('>>>>>ENZYME>>>>>',line)




        return goslim , interpro2go               

    def parse_eggnog_emapper(args,kodict):
        eggfile = os.path.join(args.v,'eggnog/nc_genomes_proteins.emapper.annotations')
        #eggfile = 'nc_genomes_proteins.emapper.annotations'
        
        COG_groups = {'A':('INFORMATION STORAGE AND PROCESSING','RNA processing and modification'),
        'B':('INFORMATION STORAGE AND PROCESSING','Chromatin Structure and dynamics'),
        'C':('METABOLISM','Energy production and conversion'),
        'D':('CELLULAR PROCESSES AND SIGNALING','Cell cycle control and mitosis'),
        'E':('METABOLISM','Amino Acid metabolis and transport'),
        'F':('METABOLISM','Nucleotide metabolism and transport'),
        'G':('METABOLISM','Carbohydrate metabolism and transport'),
        'H':('METABOLISM','Coenzyme metabolis'),
        'I':('METABOLISM','Lipid metabolism'),
        'J':('INFORMATION STORAGE AND PROCESSING','Translation'),
        'K':('INFORMATION STORAGE AND PROCESSING','Transcription'),
        'L':('INFORMATION STORAGE AND PROCESSING','Replication and repair'),
        'M':('CELLULAR PROCESSES AND SIGNALING','Cell wall/membrane/envelop biogenesis'),
        'N':('CELLULAR PROCESSES AND SIGNALING','Cell motility'),
        'O':('CELLULAR PROCESSES AND SIGNALING','Post-translational modification, protein turnover, chaperone functions'),
        'P':('METABOLISM','Inorganic ion transport and metabolism'),
        'Q':('METABOLISM','Secondary Structure'),
        'T':('CELLULAR PROCESSES AND SIGNALING','Signal Transduction'),
        'U':('CELLULAR PROCESSES AND SIGNALING','Intracellular trafficing and secretion'),
        'Y':('CELLULAR PROCESSES AND SIGNALING','Nuclear structure'),
        'Z':('CELLULAR PROCESSES AND SIGNALING','Cytoskeleton'),
        'R':('POORLY CHARACTERIZED','General Functional Prediction only'),
        'V':('CELLULAR PROCESSES AND SIGNALING','Defense mechanisms'),
        'S':('POORLY CHARACTERIZED','Function Unknown')}
        egg_protein = {}
        if not os.path.exists(eggfile):
            print('You need to generate this file first:', eggfile)
            sys.exit(1)
        else:
            with open(eggfile,'r') as infile:
                for rawline in infile:
                    if rawline[0] == '#':
                        continue
                    line = rawline.strip().split('\t')
                    protein = line[0]
                    score = line[3]
                    preferredname, EC, konum, biGG, cazy , COG_functional_category, COG_letter = None, None , None ,None , None , None, None
                    if line[5] !='':
                        preferredname = line[5]
                    
                    EC = [ec for ec in line if re.match(r'\d+\.\d+\.\d+\.',ec)]
                    if len(EC) >= 1:
                        EC = EC[0].split(',')                  
                    
                    ### Extract Mapped KO's and clean
                    konum = [ko for ko in line if 'ko:' in ko]                    
                    if len(konum) >= 1:
                        konum = konum[0].replace('ko:','').split(',')
                    
                    ### Translate KOs to defintions
                    kodefs = []
                    for ko in konum:
                        if ko in kodict:
                            kodefs.append(kodict[ko])

                    biGG = line[-1]
                    if biGG == 'NA|NA|NA':
                        continue

                    if line[-2] in COG_groups:
                        COG_functional_category = COG_groups[line[-2]]
                        COG_letter = line[-2]
                    binid = '_'.join(protein.split('_')[:2])
                    egg_protein[protein] = {'EC':EC,'biGG':biGG,'eggscore':score,'KO':konum,'KOdef':kodefs,'COG':(COG_letter,COG_functional_category)}
        return egg_protein
    
    def parse_interpro(args,goslim,interpro2go):

        interprofile = os.path.join(args.v,'eggnog/interproscan2.tsv')

        interpro_protein = {}
        if not os.path.exists(interprofile):
            print('You need to generate this file first:', interprofile)
            sys.exit(1)
        else:
            parsed_proteins = set()
            with open(interprofile,'r') as infile:
                for rawline in infile:
                    line = rawline.strip().split('\t')
                    protein = line[0]
                    if protein in parsed_proteins:
                        continue
                    parsed_proteins.add(protein)
                    iprid, iprdesc = None , None 
                    if len(line) > 11:
                        iprid = line[11]
                        iprdesc = line[12]
                    GO_low = None
                    if iprid in interpro2go:
                        GO_low = interpro2go[iprid]['GO_low']             
                    GO = None
                    GO_medium = None
                    GO_high = None
                    if 'GO:' in rawline:
                        GO = line[13].split('|')                        

                        ### Get higher level annotation for protein

                        for GOid in GO:
                            if GOid in goslim:
                                GO_medium , GO_high = goslim[GOid]['name'], goslim[GOid]['namespace']
                    interpro_protein[protein] = {'GO':GO,'IPR':(iprid,iprdesc),'GO_low':GO_low,'GO_medium':GO_medium,'GO_high':GO_high}
        return interpro_protein
    
    def write_GO_table(args,interpro_results):

        ### Establish Unique pairs of GO low + High level functions (not NONE)
        gos = set()
        OK_proteins = set()
        go_descriptions = {}
        for entry in interpro_results:
            go_pair = (interpro_results[entry]['GO_low'] , interpro_results[entry]['GO_high'])
            if not go_pair[0] == None and not go_pair[1] == None:
                for go in interpro_results[entry]['GO']:
                    gos.add(go)
                    OK_proteins.add(entry)
                    go_descriptions[go] = go_pair

        ### For each Bin count the numbeer of Unique Pairs 
        go_bin_counts = {}
        for proteinid in interpro_results:
            binid = '_'.join(proteinid.split('_')[:2])
            
            if interpro_results[proteinid]['GO'] == None:
                continue
            if not proteinid in OK_proteins:
                continue

            if binid not in go_bin_counts:
                go_bin_counts[binid] = {str(k):0 for k in gos}
                for go in interpro_results[proteinid]['GO']:
                    go_bin_counts[binid][go] += 1
            else:
                for go in interpro_results[proteinid]['GO']:
                    go_bin_counts[binid][go] += 1

        ### Write out counts for each bin\
        go_columns = [['GO','spec_function','general']]
        for go in go_descriptions:
            go_columns.append([go, go_descriptions[go][0],go_descriptions[go][1]])

        ### Get counts 
        go_counts = []
        for binid in go_bin_counts:
            row = [binid]
            for go in go_descriptions:
                row += [go_bin_counts[binid][go]]
            go_counts.append(row)
        mat1 = np.matrix(go_columns)
        mat2 = np.matrix(go_counts)
        mat2 = mat2.transpose()
        cmat = np.column_stack((mat1,mat2))
        filename = os.path.join(args.v ,'eggnog/GO_count.txt')
        #filename = 'GO_count.txt'
        np.savetxt(filename, cmat , delimiter='\t', newline='\n', fmt='%s')

    def write_KEGG_table(args,egg_results,kodict):

        kos = set()
        OK_proteins = set()
        for entry in egg_results:
            ko = egg_results[entry]['KO']
            if len(ko) != 0:
                for k in ko: 
                    kos.add(k)
                OK_proteins.add(entry)

        ### For each Bin count the numbeer of Unique Pairs 
        ko_bin_counts = {}
        for proteinid in egg_results:
            binid = '_'.join(proteinid.split('_')[:2])
            
            if len(egg_results[proteinid]['KO']) == 0:
                continue
            if not proteinid in OK_proteins:
                continue

            if binid not in ko_bin_counts:
                ko_bin_counts[binid] = {str(k):0 for k in kos}
                for ko in egg_results[proteinid]['KO']:
                    ko_bin_counts[binid][ko] += 1
            else:
                for go in egg_results[proteinid]['KO']:
                    ko_bin_counts[binid][ko] += 1

        ### Write out counts for each bin\
        ko_columns = [['KO','desc']]
        for ko in kos:
            desc = 'NA'
            if ko in kodict:
                desc = kodict[ko]
            ko_columns.append([ko,desc])
        
        ko_counts = []
        for binid in ko_bin_counts:
            row = [binid]
            for ko in kos:
                row += [ko_bin_counts[binid][ko]]
            ko_counts.append(row)

        mat1 = np.matrix(ko_columns)
        mat2 = np.matrix(ko_counts)
        mat2 = mat2.transpose()
        cmat = np.column_stack((mat1,mat2))
        filename = os.path.join(args.v ,'eggnog/KO_count.txt')
        #filname = 'eggnog/KO_count.txt'
        np.savetxt(filename, cmat , delimiter='\t', newline='\n', fmt='%s')

    def write_COG_table(args,egg_results):

        COG_groups = {'A':('INFORMATION STORAGE AND PROCESSING','RNA processing and modification'),
        'B':('INFORMATION STORAGE AND PROCESSING','Chromatin Structure and dynamics'),
        'C':('METABOLISM','Energy production and conversion'),
        'D':('CELLULAR PROCESSES AND SIGNALING','Cell cycle control and mitosis'),
        'E':('METABOLISM','Amino Acid metabolis and transport'),
        'F':('METABOLISM','Nucleotide metabolism and transport'),
        'G':('METABOLISM','Carbohydrate metabolism and transport'),
        'H':('METABOLISM','Coenzyme metabolis'),
        'I':('METABOLISM','Lipid metabolism'),
        'J':('INFORMATION STORAGE AND PROCESSING','Translation'),
        'K':('INFORMATION STORAGE AND PROCESSING','Transcription'),
        'L':('INFORMATION STORAGE AND PROCESSING','Replication and repair'),
        'M':('CELLULAR PROCESSES AND SIGNALING','Cell wall/membrane/envelop biogenesis'),
        'N':('CELLULAR PROCESSES AND SIGNALING','Cell motility'),
        'O':('CELLULAR PROCESSES AND SIGNALING','Post-translational modification, protein turnover, chaperone functions'),
        'P':('METABOLISM','Inorganic ion transport and metabolism'),
        'Q':('METABOLISM','Secondary Structure'),
        'T':('CELLULAR PROCESSES AND SIGNALING','Signal Transduction'),
        'U':('CELLULAR PROCESSES AND SIGNALING','Intracellular trafficing and secretion'),
        'Y':('CELLULAR PROCESSES AND SIGNALING','Nuclear structure'),
        'Z':('CELLULAR PROCESSES AND SIGNALING','Cytoskeleton'),
        'R':('POORLY CHARACTERIZED','General Functional Prediction only'),
        'V':('CELLULAR PROCESSES AND SIGNALING','Defense mechanisms'),
        'S':('POORLY CHARACTERIZED','Function Unknown')}

        ### For each Bin count the numbeer of Unique Pairs 
        COG_bin_counts = {}
        for proteinid in egg_results:
            binid = '_'.join(proteinid.split('_')[:2])
            
            if egg_results[proteinid]['COG'][0] is None:
                continue

            if binid not in COG_bin_counts:
                COG_bin_counts[binid] = {str(k):0 for k in COG_groups}
                cog = egg_results[proteinid]['COG'][0]
                COG_bin_counts[binid][cog] += 1
            else:
                cog = egg_results[proteinid]['COG'][0]
                COG_bin_counts[binid][cog] += 1

        ### Write out COG category counts 
        ko_columns = [['COG','High_desc','Low_desc']]
        for k in COG_groups:
            ko_columns.append([k,COG_groups[k][0],COG_groups[k][1]])
        
        ko_counts = []
        for binid in COG_bin_counts:
            row = [binid]
            for k in COG_groups:
                row += [COG_bin_counts[binid][k]]
            ko_counts.append(row)

        mat1 = np.matrix(ko_columns)
        mat2 = np.matrix(ko_counts)
        mat2 = mat2.transpose()
        cmat = np.column_stack((mat1,mat2))
        filename = os.path.join(args.v ,'eggnog/COG_count.txt')
        #filename = 'COG_count.txt'
        np.savetxt(filename, cmat , delimiter='\t', newline='\n', fmt='%s')

    def write_IPR_table(args,interpro_results):
        
        ### Establish Unique pairs of GO low + High level functions (not NONE)
        IPRs = set()
        OK_proteins = set()
        IPR_desc = {}
        for proteinid in interpro_results:
            ipr = interpro_results[proteinid]['IPR'][0]
            if ipr is None:
                continue
            IPRs.add(ipr)
            IPR_desc[ipr] = interpro_results[proteinid]['IPR'][1]
            OK_proteins.add(proteinid)
        
        IPR_bin_counts = {}
        for proteinid in interpro_results:
            binid = '_'.join(proteinid.split('_')[:2])
            
            if not proteinid in OK_proteins:
                continue
            if binid not in IPR_bin_counts:
                IPR_bin_counts[binid] = {str(k):0 for k in IPRs}
                ipr = interpro_results[proteinid]['IPR'][0]
                if ipr in IPRs:
                    IPR_bin_counts[binid][ipr] += 1 
            else:
                ipr = interpro_results[proteinid]['IPR'][0]
                if ipr in IPRs:
                    IPR_bin_counts[binid][ipr] += 1 
        
        ### Write out IPR counts 
        ipr_columns = [['IPR','IPR_desc']]
        for k in IPRs:
            ipr_columns.append([k,IPR_desc[k]])
        
        ipr_counts = []
        for binid in IPR_bin_counts:
            row = [binid]
            for k in IPRs:
                row += [IPR_bin_counts[binid][k]]
            ipr_counts.append(row)

        mat1 = np.matrix(ipr_columns)
        mat2 = np.matrix(ipr_counts)
        mat2 = mat2.transpose()
        cmat = np.column_stack((mat1,mat2))
        filename = os.path.join(args.v ,'eggnog/IPR_count.txt')
        #filename = 'IPR_count.txt'
        np.savetxt(filename, cmat , delimiter='\t', newline='\n', fmt='%s')
        
    def write_GPs(args,interpro_results):

        ### Load in Genome Property DB
        def load_GP_db():
            gp_file = '/home/projects/cpr_10006/projects/phamb/tools/genome-properties/flatfiles/genomeProperties.txt'
            GP = dict()
            AC = None
            with open(gp_file,'r') as infile:
                for line in infile:
                    if line[:2] == '--':
                        continue
                    if line[:2] == '//':
                        continue

                    line = line.strip().split('  ')
                    if len(line) > 2:
                        continue
                    entry, info = line

                    if entry == 'AC':
                        AC = info
                        GP[AC] = dict()
                        GP[AC]['TH'] = None
                        GP[AC]['TP'] = None
                        GP[AC]['DE'] = None
                        GP[AC]['current_step'] = None
                        GP[AC]['Steps'] = dict()
                        GP[AC]['AllIPR'] = None
                        GP[AC]['nrequired_steps'] = 0
                    elif entry == 'TH':
                        GP[AC]['TH'] = int(info)
                    elif entry == 'TP':
                        GP[AC]['TP'] = info
                    elif entry == 'DE':
                        GP[AC]['DE'] = info
                    elif entry == 'SN':
                        GP[AC]['Steps'][int(info)] = set()
                        GP[AC]['current_step'] = int(info)
                        GP[AC]['nrequired_steps'] += 1
                    elif entry == 'RQ':
                        current_step = GP[AC]['current_step'] 
                        if entry == 0:
                            GP[AC]['Steps'][current_step] = 'Not_Required'
                            GP[AC]['nrequired_steps'] -= 1
                    elif entry == 'EV':
                        current_step = GP[AC]['current_step'] 
                        evs = info.split(';')
                        if GP[AC]['Steps'][current_step] == 'Not_Required':
                            continue
                        else:
                            IPR = evs[0]
                            GP[AC]['Steps'][current_step].add(IPR)
            
            ### For each GP make a superset of all IPRids - if a genome does not have a single one of them, calculation for that GP can be skipped
            for gp in GP:
                iprs = set()
                for step in GP[gp]['Steps']:
                    iprs= iprs.union(GP[gp]['Steps'][step])
                GP[gp]['AllIPR'] = iprs

            return GP 
            
        def calculate_bin_GPs(GP,iprs):
            '''Evaluate if a Genome fills out any Genome Properties Partially, Complete or not at all
            Not at all = 0 
            Partially = 0.5
            Complete = 1
            '''
            gp_result = {str(k):0 for k in GP}
            for gp in GP:
                decision = 0
                all_gp_iprs = GP[gp]['AllIPR']
                if len(iprs.intersection(all_gp_iprs)) == 0:
                    continue
                Threshold = GP[gp]['TH']
                steps_completed = 0
                for step in GP[gp]['Steps']:
                    step_iprs = GP[gp]['Steps'][step]
                    if len(iprs.intersection(step_iprs))>0:
                        steps_completed += 1
                if steps_completed > Threshold:
                    decision = 0.5
                    if steps_completed >= GP[gp]['nrequired_steps']:
                        decision = 1
                gp_result[gp] = decision
            return gp_result

        ### Load Genome Property DB 
        GP = load_GP_db()

        ### Fetch all IPR entries recorded in each Genome/Bin
        IPR_bins = dict()
        for proteinid in interpro_results:
            binid = '_'.join(proteinid.split('_')[:2])
            ipr = interpro_results[proteinid]['IPR'][0]
            if ipr is None:
                continue
            if binid not in IPR_bins:
                IPR_bins[binid] = set([ipr])
            else:
                IPR_bins[binid].add(ipr)

        ### Parse all Interpros of each Genome and evaluate if any GPs are completed.
        GP_bins = dict()
        for binid in IPR_bins:
            iprs = IPR_bins[binid]
            gp_result = calculate_bin_GPs(GP,iprs)
            GP_bins[binid] = gp_result
    

        ### Finally write out matrix with the results.
        GP_columns = [['GP','GP_TP','GP_DE']]
        for gp in GP:
            GP_columns.append([gp,GP[gp]['TP'],GP[gp]['DE']])
        
        GP_counts = []
        for binid in GP_bins:
            row = [binid]
            for gp in GP:
                row += [ GP_bins[binid][gp] ]
            GP_counts.append(row)

        mat1 = np.matrix(GP_columns)
        mat2 = np.matrix(GP_counts)
        mat2 = mat2.transpose()

        ### Identify zero-rows 
        boolean = np.array([True for i in range(mat2.shape[0]) ])
        for i,row in enumerate(mat2):
            if i == 0:
                continue
            rowsum = np.array(row,dtype='float').sum()
            if rowsum == 0:
                boolean[i] = False
        mat2 = mat2[boolean]
        mat1 = mat1[boolean]
        
        ### sort 
        cmat = np.column_stack((mat1,mat2))
        filename = os.path.join(args.v ,'eggnog/GP_count.txt')
        #filename = 'GP_count.txt'
        np.savetxt(filename, cmat , delimiter='\t', newline='\n', fmt='%s')


    ### main steps
    kodict = load_KO_annotation()
    goslim,interpro2go = load_GO_mapping()

    ### Parse annotation files 
    egg_results = parse_eggnog_emapper(args,kodict)
    interpro_results = parse_interpro(args, goslim, interpro2go)

    print('Done loading annotations - starting writing')


    ### write tables
    write_GO_table(args,interpro_results)
    write_GPs(args,interpro_results)
    write_KEGG_table(args,egg_results,kodict)
    write_COG_table(args,egg_results)
    write_IPR_table(args,interpro_results)



def organise_AMRs(args):

    def CARD_RGI_enrichment_test(args):


        testdir = '/home/projects/cpr_10006/projects/phamb/HMP2/11_CARD_enrichment_analysis/'
        samplefile = os.path.join(testdir,'100samples.txt')
        samples = set()
        with open(samplefile,'r') as infile:
            for line in infile:
                line = line.strip()
                samples.add(line)

        sample_hits = []
        for sample in samples:
            rgi_file = os.path.join(testdir,'rgi',sample+'.txt')
            ngood_hits = 0
            nhits = 0
            with open(rgi_file,'r') as infile:
                infile.readline()
                for line in infile:
                    line = line.strip().split('\t')
                    nhits += 1
                    proteinid = line[1]
                    besthit_identity = float(line[9])
                    query_coverage = float(line[20])
                    if query_coverage >= 65 and besthit_identity >= 50:
                        ngood_hits += 1
            sample_hits.append(ngood_hits/2000)


    def AMRs(args):

        class make_args(object):
            def __init__(self):
                self.v = '/home/projects/cpr_10006/projects/phamb/HMP2/07_binannotation/checkv/bins_checkv_out'
    
        args = make_args()

        def calculate_stuff():

            ### Count the number of genes - sum n. genes 

            classcount = dict()
            binhits = dict()
            for proteinid in CARDhits:
                binid = '_'.join(proteinid.split('_')[:2])
                if not binid in binhits:
                    binhits[binid] = set([proteinid])
                else:
                    binhits[binid].add(proteinid)
                clas = CARDhits[proteinid][2]
                if not clas in classcount:
                    classcount[clas] = 1
                else:
                    classcount[clas] += 1

            clusterhits = dict()
            ngenes = 0
            for binid in binhits:
                cluster = binid.split('_')[1]
                if not cluster in clusterhits:
                    clusterhits[cluster] = binhits[binid]
                    ngenes += len(binhits[binid])
            print(ngenes)
            ### ~2000 
        
        def load_Resfams_db():
            ### Load DB of Resfams
            databasedirectory = '/home/projects/cpr_10006/projects/phamb/databases/'
            db = os.path.join(databasedirectory,'AMR/Resfams-full.hmm')
            AMRcutoff = dict()
            ACC = None
            with open(db,'r') as infile:
                for line in infile:
                    if line[0:3] == 'ACC':
                        ACC = line.strip().split()[1]
                        AMRcutoff[ACC] = None
                    if line[0:2] == 'TC':
                        tc = float(line.strip().split()[1])
                        AMRcutoff[ACC] = tc

            AMRannot = dict()
            resfam_annotation = os.path.join(databasedirectory,'AMR/Resfams.metadata.txt')
            with open(resfam_annotation,'r') as infile:
                for line in infile:
                        line = line.strip().split('\t')
                        ACC = line[0]
                        familyname = line[1]
                        desc = line[2]
                        mechanism = line[7]
                        AMRannot[ACC] = (familyname,mechanism,desc)
            return AMRcutoff, AMRannot
        
        def run_Resfams(args,proteinfile):
            filebase = os.path.basename(proteinfile)
            databasedirectory = '/home/projects/cpr_10006/projects/phamb/databases/'
            RESFAM_ref = os.path.join(args.v, 'AMR', filebase + '.Resfams.hmmsearch')
            if os.path.exists(RESFAM_ref):
                print('Hmmsearch file allready created:',RESFAM_ref,' delete file to recreate it')
            else:
                try:
                    db = os.path.join(databasedirectory,'AMR/Resfams-full.hmm')
                    command = [HMM_executable,
                            '--cut_ga',
                            '--cpu','24',
                            '--tblout',RESFAM_ref,
                            db,
                            proteinfile]
                    subprocess.check_call(command,stdout=subprocess.DEVNULL)
                except:
                    message = 'ERROR: Hmmsearch finished abnormally.'
                    print(message)
                    sys.exit(1)
                
                print('HMMsearch against Resfam database finished.')

        def organise_clusters(CARD_Resfam_hits,Resfamshits):

            clusterhits = dict()
            for proteinid in CARD_Resfam_hits:
                binid = '_'.join( proteinid.split('_')[:2] )
                cluster = binid.split('_')[1]
                if not cluster in clusterhits:
                    clusterhits[cluster] = dict()
                
                if not binid in clusterhits[cluster]:
                    clusterhits[cluster][binid] = set([proteinid])
                else:
                    clusterhits[cluster][binid].add(proteinid)
            return clusterhits


        def run_Blastp_CARD(args,proteinfile):
             
            blastp_executable= '/services/tools/ncbi-blast/2.8.1+/bin/blastp'
            filebase = os.path.basename(proteinfile)
            databasedirectory = '/home/projects/cpr_10006/projects/phamb/databases/'
            blastfile = os.path.join(args.v, 'AMR', filebase + '.CARDblast.m6')
            if os.path.exists(blastfile):
                print('BlastP file allready created:',blastfile,' delete file to recreate it')
            else:
                try:
                    db = os.path.join(databasedirectory,'CARD/localDB/proteindb.fsa')
                    command = [blastp_executable,
                            '-task','blastp',
                            '-evalue','0.0001',
                            '-num_threads','24',
                            '-db',db,
                            '-query',proteinfile,
                            '-out',blastfile,
                            '-outfmt', '6 std qlen slen']
                    subprocess.check_call(command)
                except:
                    message = 'ERROR: BlastP finished abnormally.'
                    print(message)
                    sys.exit(1)
            
        def parse_Blastp_hits(blastfile=None,query_coverage_cut=65,identity_cut=25,resfam = None):
            rgifile = blastfile
            CARDhits = dict()
            CARD_Resfam_hits_blast = dict()
            with open(rgifile,'r') as infile:
                for line in infile:
                    line = line.strip().split('\t')
                    proteinid = line[0]
                    card_acc = line[1]
                    identity = float(line[2])
                    query_coverage = float(line[-2])/float(line[-1])*100   
                    if query_coverage >= query_coverage_cut and identity >= identity_cut:
                        if proteinid in Resfamshits:
                            CARD_Resfam_hits_blast[proteinid] = (query_coverage,identity)
                        if not proteinid in CARDhits:
                            CARDhits[proteinid] = (query_coverage,identity)
                        else:
                            current_score = CARDhits[proteinid][1]
                            if identity >= current_score:
                                CARDhits[proteinid] = (query_coverage,identity)
            return CARDhits, CARD_Resfam_hits_blast

        def parse_Resfam_hits(hmmfile=None):
            Resfamshits  = dict()
            if not os.path.exists(hmmfile):
                print('You need to create the Resfam file? : ', hmmfile)
            else:
                with open(hmmfile,'r') as infile:
                    for line in infile:
                        if not line[0] == '#':
                            line = line.strip().split()
                            proteinid = line[0]
                            acc = line[3]
                            score = line[5]
                            if not proteinid in Resfamshits:
                                Resfamshits[proteinid] = (score,acc)
                            else:
                                current_score = Resfamshits[proteinid][0]
                                if score >= current_score:
                                    Resfamshits[proteinid] = (score,acc)   
            return Resfamshits
  
        AMRcutoff, AMRannot = load_Resfams_db()
  
        ### Write out proteins of Loose RGI-CARD hits
        rgifile = 'AMR/RGI.nucl.txt'
        rgifile_proteins_out = 'AMR/RGI.nucl.faa'
        with open(rgifile,'r') as infile,open(rgifile_proteins_out,'w') as out:
            infile.readline()
            for line in infile:
                line = line.strip().split('\t')
                proteinid = line[1].strip()
                AAsequence = line[18]
                out.write('>' + proteinid + '\n' + AAsequence + '\n')
        
        ### Run Hmmsearch towards ResFams - Also run a BlastP against the CARD database
        run_Resfams(args,'AMR/RGI.nucl.faa')
        run_Blastp_CARD(args,'AMR/RGI.nucl.faa')

        Resfamshits = parse_Resfam_hits(hmmfile = 'AMR/RGI.nucl.faa.Resfams.hmmsearch')
        blastp_CARDhits, CARD_Resfam_hits_blastp = parse_Blastp_hits(blastfile='AMR/RGI.nucl.faa.CARDblast.m6',resfam=Resfamshits)

        
        ### Parse results of RGI-CARD search
        rgifile = 'AMR/RGI.nucl.txt'
        CARDhits = dict()
        CARD_Resfam_hits = dict()
        ngoodhits = 0
        with open(rgifile,'r') as infile:
            infile.readline()
            for line in infile:
                line = line.strip().split('\t')
                proteinid = line[1].strip()
                besthit_identity = float(line[9])
                query_coverage = float(line[20])
                acc = 'AOC:'+line[10]
                classmech = line[15]
                classdesc = line[16]
                AROterm = line[8]
                print(acc,AROterm)

                if proteinid in Resfamshits:
                    CARD_Resfam_hits[proteinid] = (query_coverage,besthit_identity,classdesc,classmech,acc,AROterm)
                
                if query_coverage >= 65 and besthit_identity >= 50:
                    ngoodhits += 1
                    if not proteinid in CARDhits:
                        CARDhits[proteinid] = (query_coverage,besthit_identity,classdesc,classmech,acc,AROterm)
                    else:
                        current_score = CARDhits[proteinid][1]
                        classdesc = line[15]
                        if besthit_identity >= current_score:
                            CARDhits[proteinid] = (query_coverage,besthit_identity,classdesc,classmech,acc,AROterm)
        
        ### How many Genomes encode proteins with both a Resfam and CARD annotation?
        binhits = dict()
        for proteinid in CARD_Resfam_hits:
            binid = '_'.join( proteinid.split('_')[:2] )
            if not binid in binhits:
                binhits[binid] = set([proteinid])
            else:
                binhits[binid].add(proteinid)

        clusterhits = organise_clusters(CARD_Resfam_hits,Resfamshits)

        print(len(CARD_Resfam_hits),' Proteins annotated with both the ResFam and CARD database')
            
        ### Write out annotation from 
        #classcount = dict()
        outfile = os.path.join(args.v,'AMR/CARD.Resfams.count.txt')
        with open(outfile,'w') as out:
            out.write('cluster\tgenome\tprotein\tResfamAcc\tResfamClass\tCardARO\tCardMechanism\tCardClass\tCardSeqidentity\n')
            for cluster in clusterhits:
                for binid in clusterhits[cluster]:
                    for proteinid in clusterhits[cluster][binid]:
                        
                        ### Resfam annotation
                        acc = Resfamshits[proteinid][1]
                        Resfam_familyname, Resfam_mechanism,Resfam_description = AMRannot[acc]

                        ### CARD annotation
                        query_coverage = CARD_Resfam_hits[proteinid][0]
                        seqidentity = CARD_Resfam_hits[proteinid][1]
                        CardClass = CARD_Resfam_hits[proteinid][2]   
                        CardMechanism =  CARD_Resfam_hits[proteinid][3]    
                        CardARO =  CARD_Resfam_hits[proteinid][4]              
                        AROterm = CARD_Resfam_hits[proteinid][5]

                        if len(AROterm) >= 30:
                            AROterm = AROterm[:31]

                        outline = [cluster,binid, proteinid, acc,Resfam_mechanism, CardARO,CardMechanism,CardClass,str(seqidentity)  ]
                        out.write('\t'.join(outline)+'\n')

        ### Write out the proteins
        rgifile = 'AMR/RGI.nucl.txt'
        rgifile_proteins_out = 'AMR/RGI.nucl.Resfam.Card.fna'
    
        with open(rgifile,'r') as infile,open(rgifile_proteins_out,'w') as out:
            infile.readline()
            for line in infile:
                line = line.strip().split('\t')
                proteinid = line[1].strip()
                if proteinid in CARD_Resfam_hits:
                    DNAsequence = line[17]
                    out.write('>' + proteinid + '\n' + DNAsequence + '\n')

        fastafile ='AMR/RGI.nucl.Resfam.Card.fna'
        blastfile = 'AMR/NC.MAGS.Resfam.Card.m6'
        blast_against_MAGS(args,fastafile,blastfile)

        def blast_against_MAGS(args,fastafile,blastfile):
            
            blastp_executable= '/services/tools/ncbi-blast/2.8.1+/bin/blastn'
            
            if os.path.exists(blastfile):
                print('BlastP file allready created:',blastfile,' delete file to recreate it')
            else:
                try:
                    db = '/home/projects/cpr_10006/projects/phamb/HMP2/07_binannotation/bacteria/checkm/all_NC_genomes.fna'
                    command = [blastp_executable,
                            '-task','megablast',
                            '-perc_identity','90',
                            '-max_target_seqs','15',
                            '-evalue','0.0001',
                            '-num_threads','24',
                            '-db',db,
                            '-query',fastafile,
                            '-out',blastfile,
                            '-outfmt', '6 std qlen slen']
                    subprocess.check_call(command)
                except:
                    message = 'ERROR: BlastP finished abnormally.'
                    print(message)
                    sys.exit(1)



        













if __name__ == "__main__":
    args = parser.parse_args()

    
    nc_viral_bins_file = os.path.join(args.v,'nc_viralbins.txt')

    ### Load in ncbins 
    ncbins = set()
    with open(nc_viral_bins_file,'r') as infile:
        for line in infile:
            binid = line.strip()
            ncbins.add(binid)

    ### 
    organise_functional_annotations(args,ncbins)

    ### Antibiotic Restistance / AMRs

    organise_AMRs(args)



