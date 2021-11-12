
#!/usr/bin/python
import sys
import argparse
import vambtools as _vambtools
import os
import CAMISIM_viral_benchmark

parser = argparse.ArgumentParser(
    description="""Command-line benchmark utility.""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False)

parser.add_argument('clusterspath', help='Path to clusters.tsv')
parser.add_argument('refpath', help='Path to reference file')
parser.add_argument('directoryout', help='Path to directory of Benchmark files')
parser.add_argument('-m', dest='min_bin_size', metavar='', type=int,
                    default=2000, help='Minimum size of bins [2000]')
parser.add_argument('-s', dest='separator', help='Binsplit separator', default=None)
parser.add_argument('--disjoint', action='store_true', help='Enforce disjoint clusters')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

args = parser.parse_args()

'''
example run - assuming viral annotation and contig annotation files are available.  
python CAMISIM_viral_benchmark_cmd.py vamb_clusters.tsv pooled_gsa.txt cami_benchmark
'''

with open(args.clusterspath) as file:
    clusters = _vambtools.read_clusters(file)

with open(args.refpath) as file:
    reference = CAMISIM_viral_benchmark.Reference.from_file(file, minimum_contig_len=100)


contig_by_organism = reference.return_contig_by_organism()
binning = CAMISIM_viral_benchmark.Binning(clusters, reference, minsize=args.min_bin_size, disjoint=args.disjoint,
                            binsplit_separator=args.separator)

### Return Genome-wise performance stats for the best representative bin and write to directory out
genome_stats = binning.get_genome_performance_stats(binning.intersectionsof, contig_by_organism, clusters)
binning.write_genome_stats(genome_stats=genome_stats, directory = args.directoryout)

for rank in range(len(binning.counters)):
    binning.print_matrix_long(rank,fileout=os.path.join(args.directoryout,'RecallPrecison.tsv'))
    print("")

print('Counts straitified by organism based on; ',binning.eligible_references , 'Genomes')
for organism in binning.counters_by_organism.keys():
    binning.print_matrix_long_byorganism(organism,fileout=os.path.join(args.directoryout,'RecallPrecison.'+organism +'.tsv'))
    print("")

print('Counts based on; ',binning.eligible_references , 'Genomes')
with open(os.path.join(args.directoryout,'RecallPrecison.matrix.tsv'),'w') as out:
    for rank in range(len(binning.counters)):
        binning.print_matrix(rank,file=out)
        print("")

### Define 
viral_annotation_files = {'virusseeker':os.path.join('virus_annotations','seeker.txt'),
                    'virsorter2':os.path.join('virus_annotations','virsorter2','final-viral-score.tsv'),
                    'virfinder':os.path.join('virus_annotations','virfinder.txt'),
                    'viralverify':os.path.join('virus_annotations','viralverify','pooled_contigs.2000.tmp_result_table.csv'),
                    'deepvirfinder':os.path.join('virus_annotations','DVF','pooled_contigs.2000.tmp.fna_gt2000bp_dvfpred.txt'),
                    'voghmm':os.path.join('contig_annotations','all.hmmVOG.tbl'),
                    'micompletehmm':os.path.join('contig_annotations','all.hmmMiComplete105.tbl')
                    }

viral_annotation = CAMISIM_viral_benchmark.Viral_annotation(annotation_files=viral_annotation_files,genomes=reference)

### Benchmarks
RF_model = '../workflows/mag_annotation/dbs/RF_model.sav'
viral_benchmark = CAMISIM_viral_benchmark.Viral_benchmarking( RF_model=RF_model, contigs = viral_annotation.contigs, genomes = viral_annotation.genomes)

CAMISIM_viral_benchmark.Viral_benchmarking._write_Viraltools_benchmark_results(results_directory='viral_benchmark_bin',performance=viral_benchmark.Viraltools_benchmark_bin,bin_level=True)
CAMISIM_viral_benchmark.Viral_benchmarking._write_Viraltools_benchmark_results(results_directory='viral_benchmark_contig',performance=viral_benchmark.Viraltools_benchmark,bin_level=False)
