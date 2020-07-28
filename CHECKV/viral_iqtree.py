#!/bin/python
import argparse
import os
import sys
import pathlib
import taxopy
from Bio import SeqIO
import subprocess
import copy
import gzip

### custom modulees
import tax_annotations
import checkv_parsers

parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
''')

parser.add_argument('-c',help='VAMB clusterfile/clusters.tsv')
parser.add_argument('-a',help='Annotation summary files directory')
parser.add_argument('-v', help='checkv directory')


### Orthologous proteins for making Phylogenetic Tree
def ortholog_proteins(args):
    



        def get_terl_proteins(checkv_directory):
            '''
            VOG terminase markers are found with a simple grep-command: " grep -i terminase vog.annotations.tsv  | grep -i large  " 
            '''
            terl_proteins = set()
            terl_genomes = dict()
            
            ### Annotation ID's 

            VOG = set(['VOG00419','VOG00699','VOG00709','VOG00731','VOG00732','VOG01032','VOG01094','VOG01180','VOG01426'])
            crassterminase = 'YP_009052554.1'


            ### Annotation files to parse 
            blastp_crass = os.path.join(checkv_directory,'nc_genomes.crassphage.polterm.m6')
            voghmmfile = os.path.join(checkv_directory,'nc_genomes.VOG.hmmsearch')
            eggfile = os.path.join(checkv_directory,'eggnog/nc_genomes_proteins.emapper.annotations')
            interprofile = os.path.join(checkv_directory,'eggnog/interproscan.tsv')

            ### Crassphage markers
            with open(blastp_crass,'r') as infile:
                for line in infile:
                    line = line.strip().split('\t')
                    proteinid = line[0]
                    target = line[1]
                    binid = '_'.join(proteinid.split('_')[:2])
                    aln = int(line[3])
                    evalue = float(line[10])
                    if aln < 350:
                        continue
                    if target == crassterminase:
                        terl_proteins.add(proteinid)
                        if not binid in terl_genomes:
                            terl_genomes[binid] = set([target])
                        else:
                            terl_genomes[binid].add(target)
            
            ### VOG markers 
            with open(voghmmfile,'r') as infile:
                for line in infile:
                    if line[0] =='#':
                        continue
                    line = line.strip().split(' ')
                    line = [ x for x in line if x != '' ]
                    proteinid = line[0]
                    binid = '_'.join(proteinid.split('_')[:2])
                    vog = line[2]
                    score = float(line[5])
                    if score >= 30 and vog in VOG:
                        terl_proteins.add(proteinid)
                        if not binid in terl_genomes:
                            terl_genomes[binid] = set([vog])
                        else:
                            terl_genomes[binid].add(vog)

            ### Fluffy String Search for Interpro and Eggnog
            ### Eggnog
            with open(eggfile,'r') as infile:
                for line in infile:
                    if line[0] == '#':
                        continue
                    line = line.strip().split('\t')
                    proteinid = line[0]
                    egg = line[1]
                    binid = '_'.join(proteinid.split('_')[:2])
                    eggdescription = line[-1].lower()
                    if 'terminase' in eggdescription and 'large' in eggdescription:
                        terl_proteins.add(proteinid)
                        if not binid in terl_genomes:
                            terl_genomes[binid] = set([egg])
                        else:
                            terl_genomes[binid].add(egg)
            ### Interpro 
            with open(interprofile,'r') as infile:
                for line in infile:
                    line = line.strip().split('\t')
                    proteinid = line[0] 
                    binid = '_'.join(proteinid.split('_')[:2])
                    IPRid = line[-2]
                    iprdesc = line[-1].lower()
                    if 'terminase' in iprdesc and 'large' in iprdesc:
                        terl_proteins.add(proteinid)
                        if not binid in terl_genomes:
                            terl_genomes[binid] = set([IPRid])
                        else:
                            terl_genomes[binid].add(IPRid)
            return terl_proteins, terl_genomes

        def write_out_proteins(checkv_directory,marker_proteins,genome_taxonomy,marker,viralfamily=None):
            print(viralfamily)
            prodigalfile_AA = os.path.join(checkv_directory,'nc_genomes_proteins.faa')
            if not viralfamily is None:
                protein_outfile = os.path.join(checkv_directory,'orthologs',marker +'.'+ viralfamily + '.faa')
                written_proteins = set()
                outhandle = open(protein_outfile , 'w')
                for record in SeqIO.parse(open(prodigalfile_AA, 'r'), 'fasta'):
                    recordname = record.id
                    genomeid = '_'.join(recordname.split('_')[:2])
                    if recordname in marker_proteins and genomeid in genome_taxonomy:
                        if recordname in written_proteins:
                            continue
                        recordtaxonomy = genome_taxonomy[genomeid][0]
                        if recordtaxonomy == viralfamily:
                            written_proteins.add(recordname)
                            record.description = ''
                            SeqIO.write(record, outhandle, 'fasta')
                
            else:
                protein_outfile = os.path.join(checkv_directory,'orthologs',marker + '.faa')
                written_proteins = set()
                i = 0
                print(protein_outfile)
                outhandle = open(protein_outfile , 'w')
                for record in SeqIO.parse(open(prodigalfile_AA, 'r'), 'fasta'):
                    recordname = record.id
                    if recordname in marker_proteins:
                        if recordname in written_proteins:
                            continue
                        written_proteins.add(recordname)
                        SeqIO.write(record, outhandle, 'fasta')
                        i += 1

        def make_phylo_Tree(checkv_directory,marker,viralfamily=None):

            '''
            load the following: mafft/7.453 trimal/1.4.1 iqtree/1.6.8
             
            Prepare Phylofile 
            Run IQtree
            '''
            if not viralfamily is None:
                protein_outfile = os.path.join(checkv_directory,'orthologs',marker + '.'  + viralfamily +'.faa')
            else:
                protein_outfile = os.path.join(checkv_directory,'orthologs',marker + '.faa')

            alnfile = protein_outfile.replace('faa','aln')
            
            ### Run Mafft
            try:
                mafftexecutable = '/services/tools/mafft/7.453/bin/mafft'
                with open(alnfile,'w') as outfile:
                    command = [mafftexecutable,
                            '--auto', protein_outfile]
                    subprocess.run(command, stdout= outfile)
            except:
                message = 'ERROR: MAFFT finished abnormally.'
                print(message)
                sys.exit(1)


            ### Exchange X to - 
            # This command screws with the protein names... 
            #alnfilex = alnfile.replace('aln','xaln')
            #ps = subprocess.Popen(['cat',alnfile], stdout = subprocess.PIPE)
            #with open(alnfilex,'w') as outfile:
            #    command = ['tr',''''X'''',"'-'"]
            #    subprocess.run( command, stdin=ps.stdout, stdout=outfile)

            ### Run Trimal 
            try:
                trimalexecutable = '/services/tools/trimal/1.4.1/bin/trimal'
                phyfile = protein_outfile.replace('faa','phy')
                command = [trimalexecutable,
                        '-in', alnfile,
                        '-phylip',
                        '-out', phyfile]
                subprocess.check_call(command)
            except:
                message = 'ERROR: Trimal finished abnormally.'
                print(message)
                sys.exit(1)

            # module load iqtree/1.6.8
            # iqtree -s terl.crAss-like.phy -m TEST -nt 24 -safe -bb 1000 -alrt 1000 -redo
            # Best-fit model: VT+F+G4 chosen according to BIC
            # iqtree -s terl.phy -m LG -nt 24
            #try:
            #    iqtree_executable = '/services/tools/iqtree/1.6.8/bin/iqtree'
            #    command = [iqtree_executable,
            #            '-iqtree',
            #            '-s', phyfile,
            #            '-m', 'LG',
            #            '-nt','24']
            #    subprocess.check_call(command)
            #except:
            #    message = 'ERROR: IQtree finished abnormally.'
            #    print(message)
            #    sys.exit(1)

        def write_tree_annotation(checkv_directory,marker,terl_genomes,terl_proteins):

            tree_annotation_file = os.path.join(checkv_directory,'orthologs',marker + '.annotation.txt')


            
            with open(tree_annotation_file,'w') as out:
                out.write('proteinid\tgenomeid\tcluster\tfamily\tVOGclade\n')
                for proteinid in terl_proteins:
                    genomeid = '_'.join(proteinid.split('_')[:2])
                    family, VOGclade = genome_taxonomy[genomeid]
                    cluster = genomeid.split('_')[1]
                    out.write('{}\t{}\t{}\t{}\t{}\n'.format(proteinid,genomeid,cluster,family,VOGclade))
            
        
        

        
        ### Get all Putative Terminase (large subunit) proteins 
        checkv_directory = args.v
        terl_proteins, terl_genomes = get_terl_proteins(checkv_directory)
        
        ### Load Taxonomy
        genome_taxonomy = dict()
        with open(os.path.join(checkv_directory,'checkv_refined_taxonomy.txt') ,'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                genomeid = line[0]
                if genomeid in terl_genomes:
                    tax = line[2]
                    VOGclade = line[4]
                    family = tax.split(';')[4]
                    genome_taxonomy[genomeid] = [family,VOGclade]


        write_out_proteins(checkv_directory,terl_proteins,genome_taxonomy,'terl',viralfamily='crAss-like')
        make_phylo_Tree(checkv_directory, marker = 'terl',viralfamily='crAss-like')
        write_tree_annotation(checkv_directory,  'terl', terl_genomes, terl_proteins)


