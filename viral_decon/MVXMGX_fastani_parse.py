#!/bin/python 
import os
import sys


def parse_write_fastani(fastanifile,outfile,MVX_checkv, MGX_checkv,best_representative_dict):

    MGX_bins_hits = set()
    with open(fastanifile,'r') as infile, open(outfile,'w') as out:
        outheader = ['query','reference','ani','bifrag','totalfrag','fastani_coverage','genome_size_coverage','query_quality','ref_quality']
        out.write('\t'.join(outheader)+'\n')
        for line in infile:
            query, reference, ani, bifragment, totalfragments = line.strip().split()

            query_coverage = round(int(bifragment)/int(totalfragments),2)            
            if float(ani) >=80:
                query = os.path.basename(query).replace('.fna','')                
                reference = os.path.basename(reference)
                reference = reference.replace('.fna','')
                reference = '_'.join(reference.split('_')[:2])

                reference_cluster = reference.split('_')[1]

                if not reference_cluster in best_representative_dict:
                    continue
                if not query in MVX_checkv:
                    continue

                MGX_bins_hits.add(reference)
                ref_quality = MGX_checkv[reference][0]
                ref_size = MGX_checkv[reference][1]
                query_quality = MVX_checkv[query][0]
                query_size = MVX_checkv[query][1]
                genome_size_coverage = round(int(ref_size)/int(query_size),2)

                lineout = [query, reference, ani, bifragment, totalfragments, str(query_coverage), str(genome_size_coverage), query_quality, ref_quality] 
                out.write('\t'.join(lineout) + '\n')
    

    quality_scores = {5:'Complete',4:'High-quality',3:'Medium-quality',2:'Low-quality',1:'Not-determined'}
    nohits_file = outfile.replace('.tsv','.MGXbins.txt')
    with open(nohits_file,'w') as out:
        outheader = ['binid','checkv_quality','genome_size','contamination','hit']
        out.write('\t'.join(outheader)+'\n')
        for clusterid in best_representative_dict:
            binid, quality_score,genome_length,contaminaton = best_representative_dict[clusterid]
            if binid in MGX_bins_hits:
                hit = 'Hit'
            else:
                hit = 'noHit'
            lineout = [binid, quality_scores[quality_score],str(genome_length),str(contaminaton) , hit]
            out.write('\t'.join(lineout) + '\n')
            
                

def write_checkv_summary_nonredundant(outfile,MGX_checkv_file):
    '''
    Parses checkV file of VAMB bins and writes out the lines of the highest quality Viral genome of each Cluster.

    '''
    checkv_file = dict()
    checkv_dict = dict()
    checkv_header = None
    with open(MGX_checkv_file,'r') as infile:
        checkv_header = infile.readline()
        header = checkv_header.strip().split('\t')
        genome_copies_index = header.index('genome_copies')
        contamination_index = header.index('contamination')
        quality_index = header.index('checkv_quality')
        method_index = header.index('completeness_method')
        miuvig_quality = header.index('miuvig_quality')
        viral_genes_index = header.index('viral_genes')
        prophage_index = header.index('prophage')
        completeness_index = header.index('completeness')

        for fileline in infile:
            line = fileline.strip().split('\t')
            binid = line[0]

            quality = line[quality_index]
            genome_copies = float(line[genome_copies_index])
            contaminaton = float(line[contamination_index])

            if genome_copies >= 1.25:
                continue
                                    
            genome_length = float(line[1])
            checkv_dict[binid] = (quality,genome_length,contaminaton)

    ### Filter to the best represenative of each Cluster.

    quality_cluster_count = {'Complete':[],'High-quality':[],'Medium-quality':[],'Low-quality':[],'Not-determined':[]}
    quality_scores = {'Complete':5,'High-quality':4,'Medium-quality':3,'Low-quality':2,'Not-determined':1}
    best_representative_dict = dict()
    for binid in checkv_dict:
        quality,genome_length, contamination = checkv_dict[binid]
        quality_score = quality_scores[quality]
        clusterid = binid.split('_')[1]
        quality_cluster_count[quality].append(clusterid)  
        if not clusterid in best_representative_dict:
            best_representative_dict[clusterid] = (binid, quality_score,genome_length,contamination )
        else:
            binid, current_quality_score,current_genome_length,current_contaminaton  = best_representative_dict[clusterid]
            
            if quality_score > current_quality_score:
                best_representative_dict[clusterid] = (binid, quality_score,genome_length,contaminaton )
            elif genome_length > current_genome_length and contamination < current_contaminaton:
                best_representative_dict[clusterid] = (binid, quality_score,genome_length,contaminaton )
    
    ### Print Number of clusters in each Quality tier
    for qual in quality_cluster_count:
        counts = len(set(quality_cluster_count[qual]))
        print(qual,counts)

    quality_scores = {5:'Complete',4:'High-quality',3:'Medium-quality',2:'Low-quality',1:'Not-determined'}

    with open(outfile,'w') as out:
        outheader = ['binid','checkv_quality','genome_size','contamination','hit']
        out.write('\t'.join(outheader)+'\n')
        for clusterid in best_representative_dict:
            binid, quality_score,genome_length,contaminaton = best_representative_dict[clusterid]
            hit = 'NA'
            lineout = [binid, quality_scores[quality_score],str(genome_length),str(contaminaton) , hit]
            out.write('\t'.join(lineout) + '\n')

    #with open(outfile,'w') as out:
    #    out.write(checkv_header)
    #    for clusterid in best_representative_dict:
    #        binid, quality_score,genome_length,contaminaton = best_representative_dict[clusterid]
    #        fileline = checkv_file[binid]
    #        out.write(fileline)



def parse_write_AllVsAll(fastanifile,outfile, MGX_checkv):
    
    with open(fastanifile,'r') as infile, open(outfile,'w') as out:
        outheader = ['query','reference','qcluster','rcluster', 'ani','bifrag','totalfrag','fastani_coverage','genome_size_coverage','query_quality','ref_quality']
        out.write('\t'.join(outheader)+'\n')
        for line in infile:
            query, reference, ani, bifragment, totalfragments = line.strip().split()

            if float(ani) >=80:
                query_coverage = round(int(bifragment)/int(totalfragments),2)
                query = os.path.basename(query).replace('.fna','')

                
                reference = os.path.basename(reference).replace('.fna','')
                reference_genome = reference.replace('.fna','')
                reference = '_'.join(reference.split('_')[:2])
                
                if not reference in MGX_checkv or not query in MGX_checkv:
                    continue

                ref_quality = MGX_checkv[reference_genome][0]
                ref_size = MGX_checkv[reference_genome][1]
                query_quality = MGX_checkv[query][0]
                query_size = MGX_checkv[query][1]
                genome_size_coverage = round(int(ref_size)/int(query_size),2)

                reference_cluster = reference.split('_')[1]
                query_cluster = query.split('_')[1]
                lineout = [query, reference, query_cluster,  reference_cluster,ani, bifragment, totalfragments, str(query_coverage), str(genome_size_coverage), query_quality, ref_quality] 
                out.write('\t'.join(lineout) + '\n')



def parse_checkv_quality(checkv_quality_file,NR):
    '''
    NR variable is used to indicate 
    '''
    

    def removekey(d, key):
        r = dict(d)
        del r[key]
        return r


    checkv_dict = dict()
    with open(checkv_quality_file,'r') as infile:
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

            quality = line[quality_index]
            genome_copies = float(line[genome_copies_index])
            contaminaton = float(line[contamination_index])

            if genome_copies >= 2:
                continue
                                    
            genome_length = line[1]
            checkv_dict[binid] = (quality,genome_length,contaminaton)

    ### Filter to the best represenative of each Cluster.
    best_representative_dict = None
    quality_scores = {'Complete':5,'High-quality':4,'Medium-quality':3,'Low-quality':2,'Not-determined':1}
    if NR:
        best_representative_dict = dict()
        for binid in checkv_dict:
            quality,genome_length, contamination = checkv_dict[binid]
            quality_score = quality_scores[quality]
            clusterid = binid.split('_')[1]
            if clusterid not in best_representative_dict:
                best_representative_dict[clusterid] = (binid, quality_score,genome_length,contamination  )
            else:
                binid, current_quality_score,current_genome_length,current_contaminaton  = best_representative_dict[clusterid]
                
                if quality_score > current_quality_score:
                    best_representative_dict[clusterid] = (binid, quality_score,genome_length,contaminaton )
                elif genome_length > current_genome_length and contamination < current_contaminaton:
                    best_representative_dict[clusterid] = (binid, quality_score,genome_length,contaminaton )        
    
    return checkv_dict, best_representative_dict


def write_fastaniAllVsAll_table_extended(outfile,mgx_checkv,fastani_table_file):
    '''
    Filter Fastani hits and extend with additional query and reference information
    '''

    ### Get Quality Scores from checkV files of both MVX and MGX 
    print('Loading CheckV files')
    MGX_checkv,placeholder = parse_checkv_quality(mgx_checkv, NR=False)

    print('Parsing FastANI')
    parse_write_AllVsAll(fastani_table_file,outfile,MGX_checkv)


def write_fastani_table_extended(outfile,mgx_checkv,mvx_checkv,fastani_table_file):
    '''
    Filter Fastani hits and extend with additional query and reference information
    '''

    ### Get Quality Scores from checkV files of both MVX and MGX 
    print('Loading CheckV files')
    MVX_checkv, placeholder = parse_checkv_quality(mvx_checkv, NR=False)
    MGX_checkv,best_representative_dict = parse_checkv_quality(mgx_checkv, NR=True)

    print('Parsing FastANI')
    parse_write_fastani(fastani_table_file,outfile, MVX_checkv, MGX_checkv,best_representative_dict)

    

datasets = ['copsac','diabimu_t1d','HMP2']
for d in datasets:
    os.chdir('/home/projects/cpr_10006/projects/phamb/{}'.format(d))
    # MGX all vs all
    fastani = '06_fastani_MGXMVX/all.MGX_MGX.fastani.txt'
    fileout = '06_fastani_MGXMVX/all.MGX_MGX.fastani.ext.txt'
    mgx_checkv = '07_binannotation/checkv/VAMB_bins/quality_summary.tsv'

    write_fastaniAllVsAll_table_extended(fileout,mgx_checkv,fastani)




datasets = ['copsac','diabimu_t1d']
for d in datasets:
    os.chdir('/home/projects/cpr_10006/projects/phamb/{}'.format(d))
    mvx_checkv ='combined_assemblies/checkv.out/quality_summary.tsv'
    mgx_checkv = '07_binannotation/checkv/VAMB_bins/quality_summary.tsv'
    fileout = '06_fastani_MGXMVX/all.fastani_extended.tsv'
    fastani = '06_fastani_MGXMVX/mvx_mgs.all.fastani.txt'

    write_fastani_table_extended(fileout,mgx_checkv,mvx_checkv,fastani)



### 
datasets = ['HMP2','copsac','diabimu_t1d']
for d in datasets:
    print(d)
    os.chdir('/home/projects/cpr_10006/projects/phamb/{}'.format(d))
    mgx_checkv = '07_binannotation/checkv/VAMB_bins/quality_summary.tsv'
    fileout = '07_binannotation/checkv/VAMB_bins/NR.quality_summary.tsv'
    write_checkv_summary_nonredundant(fileout,mgx_checkv)


