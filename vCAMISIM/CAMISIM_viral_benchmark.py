
__doc__ = """Viral Benchmark script functioons"""

import collections as _collections
from collections import namedtuple as _namedtuple
from operator import add as _add
from itertools import product as _product
import numpy as _np
import sys as _sys
from math import sqrt as _sqrt
import vambtools as _vambtools
import os as _os

### Relevant for the RF-model prediction
from sklearn.metrics import precision_score, recall_score, roc_auc_score, roc_curve, confusion_matrix, f1_score, auc, matthews_corrcoef
from scipy import sparse
import joblib



class Contig:
    """An object representing a contig mapping to a subject at position start:end.
    Mapping positions use the half-open interval, like Python ranges and slices.
    Instantiate either with name, subject and mapping start/end:
        Contig('contig_41', 'subject_1', 11, 510)
    Or with only name and length
        Contig.subjectless('contig_41', 499)
    A subjectless Contig uses itself as a subject (implying it only maps to itself).
    """
    __slots__ = ['name', 'subject', 'start', 'end','len','contig_annotation']

    def __init__(self, name, subject, start, end):
        if end <= start:
            raise ValueError('Contig end must be higher than start')
        self.name = name
        self.subject = subject
        self.start = start
        self.end = end
        self.contig_annotation = {'voghmm':set(),'micompletehmm':set()}
        self.len = self.__len__()

    @classmethod
    def subjectless(cls, name, length):
        "Instantiate with only name and length"
        return cls(name, name, 0, length)
    
    def _add_annotation(self,annotation_type,annotation_tup):
        if annotation_type == 'voghmm':
            self.contig_annotation['voghmm'].add(annotation_tup[1])
        elif annotation_type == 'micompletehmm':
            self.contig_annotation['micompletehmm'].add(annotation_tup[1])
        else:
            self.contig_annotation[annotation_type] = annotation_tup

    def __repr__(self):
        return 'Contig({}, subject={}, {}:{})'.format(self.name, self.subject, self.start, self.end)

    def __len__(self):
        return self.end - self.start

class Genome:
    """A set of contigs known to come from the same organism.
    The breadth must be updated after adding/removing contigs with self.update_breadth(),
    before it is accurate.
    >>> genome = Genome('E. coli')
    >>> genome.add(contig)
    >>> genome.update_breadth()
    """
    __slots__ = ['name', 'breadth', 'contigs','genome_annotation']

    def __init__(self, name):
        self.name = name
        self.contigs = set()
        self.breadth = 0
        self.genome_annotation = dict()

    def add(self, contig):
        self.contigs.add(contig)

    def remove(self, contig):
        self.contigs.remove(contig)

    def discard(self, contig):
        self.contigs.discard(contig)

    @property
    def ncontigs(self):
        return len(self.contigs)

    @staticmethod
    def getbreadth(contigs):
        "This calculates the total number of bases covered at least 1x in ANY Genome."
        bysubject = _collections.defaultdict(list)
        for contig in contigs:
            bysubject[contig.subject].append(contig)

        breadth = 0
        for contiglist in bysubject.values():
            contiglist.sort(key=lambda contig: contig.start)
            rightmost_end = float('-inf')

            # We only summarise for positions covered by 1 bp (no redundant counting)
            for contig in contiglist:
                breadth += max(contig.end, rightmost_end) - max(contig.start, rightmost_end)
                rightmost_end = max(contig.end, rightmost_end)

        return breadth
    
    def update_breadth(self):
        "Updates the breadth of the genome"
        self.breadth = self.getbreadth(self.contigs)

    def __repr__(self):
        return 'Genome({}, ncontigs={}, breadth={})'.format(self.name, self.ncontigs, self.breadth)

class Reference:
    """A set of Genomes known to represent the ground truth for binning.
    Instantiate with any iterable of Genomes.
    >>> print(my_genomes)
    [Genome('E. coli'), ncontigs=95, breadth=5012521),
     Genome('Y. pestis'), ncontigs=5, breadth=46588721)]
    >>> Reference(my_genomes)
    Reference(ngenomes=2, ncontigs=100)
    Properties:
    self.genomes: {genome_name: genome} dict
    self.contigs: {contig_name: contig} dict
    self.genomeof: {contig: genome} dict
    self.breadth: Total length of all genomes
    self.ngenomes
    self.ncontigs
    """

    # Instantiate with any iterable of Genomes
    def __init__(self, genomes, taxmaps=list()):
        self.genomes = dict() # genome_name : genome dict
        self.contigs = dict() # contig_name : contig dict
        self.genomeof = dict() # contig : genome dict

        # This is a list of dicts: The first one maps genomename to name of next taxonomic level
        # The second maps name of second level to name of third level etc.
        self.taxmaps = taxmaps

        # Load genomes into list in case it's a one-time iterator
        genomes_backup = list(genomes) if iter(genomes) is genomes else genomes

        # Check that there are no genomes with same name
        if len({genome.name for genome in genomes_backup}) != len(genomes_backup):
            raise ValueError('Multiple genomes with same name not allowed in Reference.')

        for genome in genomes_backup:
            self.add(genome)

        self.breadth = sum(genome.breadth for genome in genomes_backup)

    def load_tax_file(self, line_iterator, comment='#'):
        """Load in a file with N+1 columns, the first being genomename, the next being
        the equivalent taxonomic annotation at different ranks
        Replaces the Reference's taxmaps list."""
        taxmaps = list()
        isempty = True

        for line in line_iterator:
            if line.startswith(comment):
                continue

            genomename, *clades = line[:-1].split('\t')

            if isempty:
                if not clades:
                    raise ValueError('Must have at least two columns')

                for i in clades:
                    taxmaps.append(dict())
                isempty = False

            if genomename in taxmaps[0]:
                raise KeyError("Genome name {} present more than once in taxfile".format(genomename))

            previousrank = genomename
            for nextrank, rankdict in zip(clades, taxmaps):
                existing = rankdict.get(previousrank, nextrank)
                if existing != nextrank:
                    raise KeyError("Rank {} mapped to both {} and {}".format(previousrank, existing, nextrank))

                rankdict[previousrank] = nextrank
                previousrank = nextrank

        self.taxmaps = taxmaps

    @property
    def ngenomes(self):
        return len(self.genomes)

    @property
    def ncontigs(self):
        return len(self.contigs)

    def __repr__(self):
        ranks = len(self.taxmaps) + 1
        return 'Reference(ngenomes={}, ncontigs={}, ranks={})'.format(self.ngenomes, self.ncontigs, ranks)

    @staticmethod
    def _parse_subject_line(line):
        "Returns contig, genome_name from a reference file line with subjects"
        contig_name, genome_name, taxid, contig_id,nreads, start, end = line[:-1].split('\t')
        start = int(start)
        end = int(end) + 1 # semi-open interval used in internals, like range()
        contig = Contig(contig_name, genome_name, start, end)
        return contig, genome_name

    @staticmethod
    def _parse_subjectless_line(line):
        "Returns contig, genome_name from a reference file line without subjects"
        contig_name, genome_name, length = line[:-1].split('\t')
        length = int(length)
        contig = Contig.subjectless(contig_name, length)
        return contig, genome_name

    @classmethod
    def _parse_file(cls, filehandle,subjectless=False, minimum_contig_len = 2000):
        "Returns a list of genomes from a reference file"
        function = cls._parse_subjectless_line if subjectless else cls._parse_subject_line

        genomes = dict()
        for line in filehandle:
            # Skip comments
            if line.startswith('#'):
                continue

            contig, genome_name = function(line)

            genome = genomes.get(genome_name)
            if genome is None:
                genome = Genome(genome_name)
                genomes[genome_name] = genome
            
            ### Add only contigs above a specified size (for VAMB contigs > 2000 bp)
            contig_len = contig.end - contig.start
            if contig_len >= minimum_contig_len:
                genome.add(contig)

        # Update all genomes
        genomes = list(genomes.values())
        for genome in genomes:
            genome.update_breadth()

        return genomes
    

    @classmethod
    def from_file(cls, filehandle, subjectless=False, minimum_contig_len=2000):
        """Instantiate a Reference from an open filehandle.
        "subjectless" refers to the style of reference file: If true, assumes columns are
        [contig_name, genome_name, contig_length]. If false, assume
        [contig_name, genome_name, subject_name, mapping_start, mapping_end]
        >>> with open('my_reference.tsv') as filehandle:
            Reference.from_file(filehandle)
        Reference(ngenomes=2, ncontigs=100)
        """

        genomes = cls._parse_file(filehandle, subjectless=subjectless, minimum_contig_len=minimum_contig_len)
        return cls(genomes)

    @classmethod
    def from_clusters(cls, clusters, fastadict, minimum_contig_len=2000):
        """Instantiate Genomes from VAMB clusters.
        Genome structures are obtained through VAMB clusters
        Reference(ngenomes=2, ncontigs=100)
        """

        genomes = cls._parse_clusters(clusters,fastadict,minimum_contig_len=minimum_contig_len)
        return cls(genomes)

    @classmethod
    def _parse_clusters(cls, clusters, fastadict, minimum_contig_len = 2000):
        "Returns a list of genomes from a reference file"

        genomes = dict()
        for genome_name in clusters:
            for contig in clusters[genome_name]:
                genome = genomes.get(genome_name)
                if genome is None:
                    genome = Genome(genome_name)
                    genomes[genome_name] = genome
                contig_len = fastadict[contig].__len__()
                contig = Contig(contig, genome_name, 1, contig_len)
                ### Add only contigs above a specified size (for VAMB contigs > 2000 bp)
                if contig_len >= minimum_contig_len:
                    genome.add(contig)

        # Update all genomes
        genomes = list(genomes.values())

        return genomes

    def return_contig_by_organism(self):
        '''...'''
        organism_key = {'viral':'MGV','bacteria':'GCF','plasmid':'NZ'}
        contig_by_organism = dict()
        for contig in self.genomeof:
            subject_name = contig.subject
            subject_organism = None
            for key, textkey in organism_key.items():
                if textkey in subject_name:
                    subject_organism = key
                    contig_by_organism[contig.name] = {'organism':key,'size':contig.len}
        return contig_by_organism
    
    def add(self, genome):
        "Adds a genome to this Reference. If already present, do nothing."
        if genome.name not in self.genomes:
            self.genomes[genome.name] = genome
            for contig in genome.contigs:
                if contig.name in self.contigs:
                    raise KeyError("Contig name '{}' multiple times in Reference.".format(contig.name))

                self.contigs[contig.name] = contig
                self.genomeof[contig] = genome

    def remove(self, genome):
        "Removes a genome from this Reference, raising an error if it is not present."
        del self.genomes[genome.name]

        for contig in genome.contigs:
            del self.contigs[contig.name]
            del self.genomeof[contig]

    def discard(self, genome):
        "Remove a genome if it is present, else do nothing."
        if genome.name in self.genomes:
            self.remove(genome)

    def filter_genomes(self):
        '''Remove genomes with 0 breadth i.e.not covered in dataset'''
        for genome in self.genomes.values():
            if genome.breadth == 0:
                self.remove(genome)


class Binning:
    """The result of a set of clusters applied to a Reference.
    >>> ref
    (Reference(ngenomes=2, ncontigs=5)
    >>> b = Binning({'bin1': {contig1, contig2}, 'bin2': {contig3, contig4}}, ref)
    Binning(4/5 contigs, ReferenceID=0x7fe908180be0)
    >>> b[0.5, 0.9] # num. genomes 0.5 recall, 0.9 precision
    1
    Init arguments:
    ----------- Required ---------
    contigsof:     Dict of clusters, each sequence present in the Reference
    reference:     Associated Reference object
    ----------- Optional ---------
    recalls:       Iterable of minimum recall thresholds
    precisions:    Iterable of minimum precision thresholds
    checkpresence: Whether to raise an error if a sequence if not present in Reference
    disjoint:      Whether to raise an error if a sequence is in multiple bins
    binsplit_separator: Split bins according to prefix before this separator in seq name
    minsize:       Minimum sum of sequence lengths in a bin to not be ignored
    mincontigs:    Minimum number of sequence in a bin to not be ignored
    Properties:
    self.reference:       Reference object of this benchmark
    self.recalls:         Sorted tuple of recall thresholds
    self.precisions:      Sorted tuple of precision thresholds
    self.nbins:           Number of bins
    self.ncontigs:        Number of binned contigs
    self.contigsof:       {bin_name: {contig set}}
    self.binof:           {contig: bin_name(s)}, val is str or set
    self.breadthof:       {bin_name: breadth}
    self.intersectionsof: {genome: {bin:_name: intersection}}
    self.breadth:         Total breadth of all bins
    self.counters:        List of (rec, prec) Counters of genomes for each taxonomic rank
    """
    _DEFAULTRECALLS = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
    _DEFAULTPRECISIONS = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

    @property
    def nbins(self):
        return len(self.contigsof)

    @property
    def ncontigs(self):
        return len(self.binof)

    def __init__(self, contigsof, reference, recalls=_DEFAULTRECALLS,
              precisions=_DEFAULTPRECISIONS, checkpresence=True, disjoint=True,
              binsplit_separator=None, minsize=None, mincontigs=1):
        # See class docstring for explanation of arguments

        self.precisions = tuple(sorted(precisions))
        self.recalls = tuple(sorted(recalls))
        self.reference = reference
        self.eligible_references = 0
        

        self.contigsof = dict() # bin_name: {contigs} dict
        self.binof = dict() # contig: bin_name or {bin_names} dict
        self.breadthof = dict() # bin_name: int dict
        self._parse_bins(contigsof, checkpresence, disjoint, binsplit_separator, minsize, mincontigs)
        self.breadth = sum(self.breadthof.values())

        # intersectionsof[genome] = {genome: {binname: tp, binname: tp ... }}
        # for all bins with nonzero true positives
        intersectionsof = dict()
        for genome in reference.genomes.values():
                
            if genome.breadth == 0:
                continue
            self.eligible_references += 1

            intersectionsof[genome] = dict()
            for bin_name, intersection in self._iter_intersections(genome):
                intersectionsof[genome][bin_name] = intersection
        self.intersectionsof = intersectionsof

        # Set counts
        self.counters = self._getcounts()
        self.counters_by_organism = self._getcounts_by_organism()

    def _iter_intersections(self, genome):
        """Given a genome, return a generator of (bin_name, intersection) for
        all binning bins with a nonzero recall and precision.
        """
        # Get set of all binning bin names with contigs from that genome
        bin_names = set()
        for contig in genome.contigs:
            bin_name = self.binof.get(contig)
            if bin_name is None:
                continue
            elif isinstance(bin_name, set):
                bin_names.update(bin_name)
            else:
                bin_names.add(bin_name)

        for bin_name in bin_names:
            intersecting_contigs = genome.contigs.intersection(self.contigsof[bin_name])
            intersection = Genome.getbreadth(intersecting_contigs)
            yield bin_name, intersection

    def confusion_matrix(self, genome, bin_name):
        "Given a genome and a binname, returns TP, TN, FP, FN"

        true_positives = self.intersectionsof[genome].get(bin_name, 0)
        false_positives = self.breadthof[bin_name] - true_positives
        false_negatives = genome.breadth - true_positives
        true_negatives = self.reference.breadth - false_negatives - false_positives + true_positives

        return true_positives, true_negatives, false_positives, false_negatives

    def mcc(self, genome, bin_name):
        "Calculate Matthew's correlation coefficient between a genome and a bin."

        tp, tn, fp, fn = self.confusion_matrix(genome, bin_name)
        mcc_num = tp * tn - fp * fn
        mcc_den = (tp + fp) * (tp + fn)
        mcc_den *= (tn + fp) * (tn + fn)
        return 0 if mcc_den == 0 else mcc_num / _sqrt(mcc_den)

    def f1(self, genome, bin_name):
        "Calculate F1 score between genome and a bin"

        tp, tn, fp, fn = self.confusion_matrix(genome, bin_name)
        return 2*tp / (2*tp + fp + fn)

    def _getseen(self, recprecof):
        """Make a {genome: isseen} dict, where isseen is a boolean vector
        (implemented as an integer), 1 if a genome is seen at that recall, prec level,
        0 otherwise
        """
        isseen = dict()
        for genome, _dict in recprecof.items():
            seen = 0
            for binname, (recall, precision) in _dict.items():
                for i, (min_recall, min_precision) in enumerate(_product(self.recalls, self.precisions)):
                    if recall < min_recall:
                        break

                    if precision >= min_precision:
                        seen |= 1 << i
            isseen[genome] = seen
        return isseen

    def _accumulate(self, seen, counts):
        "Given a 'seen' dict, make a dict of counts at each threshold level"
        nsums = (len(self.recalls) * len(self.precisions))
        sums = [0] * nsums
        for v in seen.values():
            for i in range(nsums):
                sums[i] += (v >> i) & 1 == 1

        for i, (recall, precision) in enumerate(_product(self.recalls, self.precisions)):
            counts[(recall, precision)] = sums[i]

    def _get_prec_rec_dict(self):
        recprecof = _collections.defaultdict(dict)
        for genome, intersectiondict in self.intersectionsof.items():
            for binname in intersectiondict:
                tp, tn, fp, fn = self.confusion_matrix(genome, binname)
                recall = tp / (tp + fn)
                precision = tp / (tp + fp)
                recprecof[genome.name][binname] = (recall, precision)

        return recprecof

    def _getcounts(self):
        # One count per rank (+1 for inclusive "genome" rank)
        counts = [_collections.Counter() for i in range(len(self.reference.taxmaps) + 1)]
        recprecof = self._get_prec_rec_dict()
        seen = self._getseen(recprecof)
        # Calculate counts for each taxonomic level
        for counter, taxmap in zip(counts, self.reference.taxmaps):
            self._accumulate(seen, counter)
            newseen = dict()
            for clade, v in seen.items():
                newclade = taxmap[clade]
                newseen[newclade] = newseen.get(newclade, 0) | v
            seen = newseen

        self._accumulate(seen, counts[-1])

        return counts

    
    def _getcounts_by_organism(self):
        '''Stratify the counts by organism.'''
        organism_key = {'viral':'MGV','bacteria':'GCF','plasmid':'NZ'}
        organism_counts = dict()
        for organism in organism_key:
            organism_counts[organism] = [_collections.Counter() for i in range(len(self.reference.taxmaps) + 1)]
        # One count per rank (+1 for inclusive "genome" rank)
        counts = [_collections.Counter() for i in range(len(self.reference.taxmaps) + 1)]
        recprecof = self._get_prec_rec_dict()
        seen = self._getseen(recprecof)
        
        # Calculate counts for each taxonomic level
        for counter, taxmap in zip(counts, self.reference.taxmaps):
            self._accumulate(seen, counter)
            newseen = dict()
            for clade, v in seen.items():
                newclade = taxmap[clade]
                newseen[newclade] = newseen.get(newclade, 0) | v
            seen = newseen

        for organism, key in organism_key.items():
            seen_organism = {k:seen[k] for k in seen if key in k}
            self._accumulate(seen_organism,organism_counts[organism][-1])
        return organism_counts

    def _parse_bins(self, contigsof, checkpresence, disjoint, binsplit_separator, minsize, mincontigs):
        "Fills self.binof, self.contigsof and self.breadthof during instantiation"

        if binsplit_separator is not None:
            contigsof = _vambtools.binsplit(contigsof, binsplit_separator)

        if minsize is not None or mincontigs is not None:
            minsize = 1 if minsize is None else minsize
            mincontigs = 1 if mincontigs is None else mincontigs
            contigsof = filter_clusters(contigsof, self.reference, minsize, mincontigs, checkpresence=checkpresence)

        for bin_name, contig_names in contigsof.items():
            contigset = set()
            # This stores each contig by their true genome name.
            contigsof_genome = _collections.defaultdict(list)

            for contig_name in contig_names:
                contig = self.reference.contigs.get(contig_name)

                # Check that the contig is in the reference
                if contig is None:
                    if checkpresence:
                        raise KeyError('Contig {} not in reference.'.format(contig_name))
                    else:
                        continue

                # Check that contig is only present one time in input
                existing = self.binof.get(contig)
                if existing is None:
                    self.binof[contig] = bin_name
                else:
                    if disjoint:
                        raise KeyError('Contig {} found in multiple bins'.format(contig_name))
                    elif isinstance(existing, str):
                        self.binof[contig] = {existing, bin_name}
                    else:
                        self.binof[contig].add(bin_name)

                contigset.add(contig)
                genome = self.reference.genomeof[self.reference.contigs[contig_name]]
                contigsof_genome[genome.name].append(contig)

            self.contigsof[bin_name] = contigset

            # Now calculate breadth of bin. This is the sum of the number of covered
            # base pairs for all the genomes present in the bin.
            breadth = 0
            for contigs in contigsof_genome.values():
                breadth += Genome.getbreadth(contigs)
            self.breadthof[bin_name] = breadth

    
    def get_genome_performance_stats(self,intersectionsof,contig_by_organism, clusters):
        '''Inspect Intersection of VAMB bin and Simulated CAMI genome'''
        genome_stats = dict()
        for genome_reference, vambbins in intersectionsof.items():
            
            if len(vambbins) == 0:
                continue

            ### sort to identify most representative bin
            vambbins_sorted = {k: v for k, v in sorted(vambbins.items(), key=lambda item: item[1])}
            best_representative_bin = list(vambbins_sorted.keys())[-1]

            organism_bp = {'viral':0,'bacteria':0,'plasmid':0} # Sum this to get the total binsize 
            for contig in clusters[best_representative_bin]:
                if contig in contig_by_organism:
                    org, contiglen = contig_by_organism[contig]['organism'], contig_by_organism[contig]['size']
                    organism_bp[org] += contiglen
            
            ### How many basepairs aligns specifically ot the reference genome
            binsize_total = sum(organism_bp.values())
            
            ### Calculate Precision Recall scores of Reference genomes 
            tp, tn, fp, fn = self.confusion_matrix(genome_reference, best_representative_bin)
            mcc_score = self.mcc(genome_reference, best_representative_bin)
            f1_score = 2*tp / (2*tp + fp + fn)
            recall = tp / (tp + fn)
            precision = tp / (tp + fp)
            genome_stats[genome_reference.name] = {'mcc':mcc_score,
            'reference_breadth':genome_reference.breadth,
            'f1':f1_score,
            'recall':recall,
            'precision':precision,
            'best_bin':best_representative_bin,
            'bin_contig_distribution':organism_bp,
            'binsize':binsize_total,
            'tp':tp, 
            'fp':fp,
            'tn':tn,
            'fn':fn}
        return genome_stats

    def write_genome_stats(self,genome_stats, directory):
        '''Write Genome recovery stats for each reference genome'''
        if not _os.path.exists(directory):
            _os.makedirs(directory)

        with open(_os.path.join(directory,'GenomeRecallStats.tsv'),'w') as out:
            header = ['Genome','genome_size','bin','binsize','tp','fp','tn','fn','mcc','f1','recall','precision','viralbp','bacteriabp','plasmidbp']
            out.write('\t'.join(header)+'\n')
            for genome in genome_stats:
                g = genome_stats[genome]
                bin_distribution = [ genome_stats[genome]['bin_contig_distribution'][k] for k in genome_stats[genome]['bin_contig_distribution']]
                row = [genome,g['reference_breadth'],g['best_bin'],g['binsize'],g['tp'],g['fp'],g['tn'],g['fn'],g['mcc'],g['f1'],g['recall'],g['precision']] + bin_distribution
                out.write('\t'.join([str(i) for i in row])+'\n')

    @classmethod
    def from_file(cls, filehandle, reference, recalls=_DEFAULTRECALLS,
                  precisions=_DEFAULTPRECISIONS, checkpresence=True, disjoint=True,
                  binsplit_separator=None, minsize=None, mincontigs=None):
        contigsof = dict()
        for line in filehandle:
            if line.startswith('#'):
                continue

            line = line.rstrip()
            bin_name, tab, contig_name = line.partition('\t')

            if bin_name not in contigsof:
                contigsof[bin_name] = [contig_name]
            else:
                contigsof[bin_name].append(contig_name)

        return cls(contigsof, reference, recalls, precisions, checkpresence, disjoint,
                   binsplit_separator, minsize, mincontigs)

    def print_matrix(self, rank, file=_sys.stdout):
        """Prints the recall/precision number of bins to STDOUT."""

        if rank >= len(self.counters):
            raise IndexError("Taxonomic rank out of range")

        print('\tRecall', file=file)
        print('Prec.', '\t'.join([str(r) for r in self.recalls]), sep='\t', file=file)

        for min_precision in self.precisions:
            row = [self.counters[rank][(min_recall, min_precision)] for min_recall in self.recalls]
            print(min_precision, '\t'.join([str(i) for i in row]), sep='\t', file=file)

    def print_matrix_long(self, rank, fileout=None):
        """Prints the recall/precision number of bins to STDOUT."""

        if rank >= len(self.counters):
            raise IndexError("Taxonomic rank out of range")

        with open(fileout,'w') as out:
            out.write('Precision\tRecall\tnGenomes\n')
            for min_precision in self.precisions:
                for min_recall in self.recalls:
                    row = [min_precision,min_recall,self.counters[rank][(min_recall, min_precision)]]
                    out.write('\t'.join([str(i) for i in row])+'\n')

    def print_matrix_long_byorganism(self, organism, fileout=None):
        """Prints the recall/precision number of bins to STDOUT."""

        with open(fileout,'w') as out:
            out.write('Precision\tRecall\tnGenomes\torganism\n')
            for min_precision in self.precisions:
                for min_recall in self.recalls:
                    row = [min_precision,min_recall,self.counters_by_organism[organism][0][(min_recall, min_precision)],organism]
                    out.write('\t'.join([str(i) for i in row])+'\n')


    def __repr__(self):
        fields = (self.ncontigs, self.reference.ncontigs, hex(id(self.reference)))
        return 'Binning({}/{} contigs, ReferenceID={})'.format(*fields)

    def summary(self, precision=0.9, recalls=None):
        if recalls is None:
            recalls = self.recalls
        return [[counter[(recall, precision)] for recall in recalls] for counter in self.counters]



class Viral_annotation:
    '''A module for parsing all kinds of annotation pr. contig for a Reference Genome or VAMB bin'''
    # Instantiate with any iterable of Genomes
    def __init__(self, genomes, annotation_files):
        '''genomes is either a Reference Genome or VAMB bin'''

        self.contigs = genomes.contigs      # Make shallow copy of contigs 
        self.genomes = genomes.genomes      # Make shallow copy of genomes

        ### Check files
        for f in annotation_files.values():
            if not _os.path.exists(f):                
                raise KeyError("One of the annotation files does not exist: {}".format(f))

        ### Parse all annotation files, add information for each contig
        for filetype, file in annotation_files.items():
            print('Parsing {}'.format(filetype))
            self._parse_viralannotation_file(filetype.lower(),file)

        ### Summarise annotation for each Genome/ aggregated contig annotation
        self._summarise_genome_annotation()

    def _parse_viralannotation_file(self,filetype,file):
        '''  '''
        function_dict = {'virusseeker':viral_annotation_tool_parsers._parse_virusseker_row,
                         'virsorter2':viral_annotation_tool_parsers._parse_vs2_row,
                         'virfinder':viral_annotation_tool_parsers._parse_virfinder_row,
                         'viralverify':viral_annotation_tool_parsers._parse_viralverify_row,
                         'deepvirfinder':viral_annotation_tool_parsers._parse_dvf_row,
                         'vibrant':viral_annotation_tool_parsers._parse_vibrant_row,
                         'voghmm':viral_annotation_tool_parsers._parse_hmm_row,
                         'micompletehmm':viral_annotation_tool_parsers._parse_hmm_row,
                         'checkv':viral_annotation_tool_parsers._parse_checkv
                         }
        with open(file,'r') as filehandle:
            header = filehandle.readline()
            if not filetype in function_dict:
                raise ValueError('Somethings wrong with the name in the File type dict')
            parse_function = function_dict[filetype]

            # We wanna add the extra annotation for each Contig but the format varies 
            for line in filehandle:
                
                if line[0] == '#':
                    continue
                annotation_tuple = parse_function(line)

                if not annotation_tuple is None:
                    contig = self.contigs[annotation_tuple[0]]
                    contig_len = contig.__len__()
                    annotation_tuple = annotation_tuple + (contig_len,)
                    contig._add_annotation(filetype,annotation_tuple)
    
    def _summarise_genome_annotation(self):
        '''Annotation is summarised over contigs for each genome'''

        for genome in self.genomes.values():

            # Establish a dict with all contigs annotations for each Annotation type i.e. Virsorter2, Viralseeker etc.
            genome_annotation_aggregated = self._aggregate_genome_annotation(genome)

            # Summarise Genome annotation for each Annotation type
            genome_annotation_aggregated_summary = dict()
            for annotation_type in genome_annotation_aggregated:
                if annotation_type == 'voghmm':
                    # The Number of Distinct VOGs / number of contigs in Genome/bin 
                    VOGscaled =  round( len( set().union(*genome_annotation_aggregated[annotation_type]) ) / genome.ncontigs , 2)
                    genome_annotation_aggregated_summary[annotation_type] = VOGscaled
                elif annotation_type == 'micompletehmm':
                    # The number of Distinct Bacterial Hallmarks in Genome/bin 
                    bacterial_hallmarks = len( set().union(*genome_annotation_aggregated[annotation_type]) )
                    genome_annotation_aggregated_summary[annotation_type] = bacterial_hallmarks
                else:
                    # Calculate aggregate prediction scores over contigs in Genome/bin 
                    contig_scores = [item[1] for item in genome_annotation_aggregated[annotation_type]]
                    contig_lengths = [item[-1] for item in genome_annotation_aggregated[annotation_type]]
                    annotated_contigs = len(contig_scores) # How many contigs actually got annotated
                    majority_vote, mean_score, weighted_mean_score, median_score = self._majority_vote(contig_scores), round(_np.average(contig_scores),2), round(_np.average(contig_scores,weights=contig_lengths),2), round(_np.median(contig_scores),2)
                    ViralResults = {'majority_vote':majority_vote,'mean_score':mean_score,'weighted_mean_score':weighted_mean_score,'median_score':median_score}
                    genome_annotation_aggregated_summary[annotation_type] = ViralResults
            
            genome.genome_annotation = genome_annotation_aggregated_summary

    def _aggregate_genome_annotation(self, genome):
        ''' Return Contig annotation for each Annotation type '''
        genome_annotation_aggregated = dict()
        for contig in genome.contigs:
            contig_name = contig.name

            contig_annotation = self.contigs[contig_name].contig_annotation

            for annotation_type in contig_annotation:
                if not annotation_type in genome_annotation_aggregated:
                    genome_annotation_aggregated[annotation_type] = []
                genome_annotation_aggregated[annotation_type] += [contig_annotation[annotation_type]]
        return genome_annotation_aggregated

    @staticmethod
    def _majority_vote(numberlist, cutoff = 0.75):
        boolean_list = _np.array(numberlist) > cutoff
        boolean = _collections.Counter(boolean_list).most_common(1)[0][0]
        return boolean


class Viral_benchmarking():
    def __init__(self,RF_model,contigs,genomes):
        self.RF_model = RF_model
        self.Viraltools_benchmark = None
        self.Viraltools_benchmark_bin = dict()
        self.genometype_lookup = dict()
        self.RF_predictions = None


        # Would be way better with a Genome metadata file
        for genome in genomes.values():
            self.genometype_lookup[genome.name] = 'bacterial' if 'GCF' in genome.name else 'viral'
        
        ### Benchmark of viral predictors on Contig-level
        self._benchmark_Viraltools(contigs,genomes)

        ### Prepare and run PHAMB
        genome_order, sparse_df = self._return_RF_prediction_dataframe(genomes)
        RF_predictions = self._run_RF_model(self.RF_model,self.genometype_lookup,genome_order,sparse_df)
        self.RF_predictions  = RF_predictions

        ### Aggregate Bin-level scores for Virus annotation tools
        metrics = ['mean_score','weighted_mean_score','median_score']
        for summary_metric in metrics:
            prediction_data_frames = self._return_viraltool_prediction_tables_binlevel(genomes,summary_metric=str(summary_metric))

            ### Determine performance
            self._benchmark_Viraltools_binlevel(prediction_data_frames,summary_metric)
        
        ### RF model
        self._benchmark_Viraltools_binlevel(RF_predictions,'RF')
        

    
    def _benchmark_Viraltools(self,contigs,genomes):
        '''Summarise the number of correctly predicted base-pairs for each Genome Type : Viral or Bacterial'''
        
        # Parse all contigs , each Viral Predictor, return the number of Predicted Bacterial and Viral base-pairs. 
        # Compare the predicted base-pairs with the ground truth (summed breadth of Bacterial and Viral genomes respectively.) 
        tools = ['virusseeker','virfinder','virsorter2','deepvirfinder','viralverify']
        viral_annotation_tools_cutoffs = {'virfinder':0.75,'virusseeker':0.75,'virsorter2':0.75,'deepvirfinder':0.75,'viralverify':5}
        
        # Breadth counts
        #breadth = {'viral':0,'bacterial':0}
        contig_types = {'viral':set(),'bacterial':set()}
        contig_type_lookup = {}
        for contig in contigs.values():
            genome_type = 'viral' if 'MGV' in contig.subject else 'bacterial'
            contig_types['viral'].add(contig) if genome_type =='viral' else contig_types['bacterial'].add(contig)
            contig_type_lookup[contig.name] = genome_type

        performance = dict()
        for tool in tools:
            predictions_table = self._return_viraltool_prediction_tables(tool,contigs,contig_type_lookup)

            # Calculate performance based on prediction scores
            auc_curve,prerec_curve, general_scores = self._calculate_TPR_FPR_AUC(predictions_table)
            performance[tool] = {'AUC':auc_curve,'PreRec':prerec_curve,'general_scores':general_scores}
            # Calculate performance based on covered BPs according to predefined thresholds 

        # Save results
        self.Viraltools_benchmark = performance


    def _benchmark_Viraltools_binlevel(self,prediction_data_frames,summary_metric):
        '''ON GENOME LEVEL: Parses a list of data-frames with predicted labels for Genomes'''
        performance = dict()
        # Calculate performance on Genome level
        if summary_metric == 'RF':
            prediction_table = prediction_data_frames
            prediction_tool = prediction_table[0][0]
            auc_curve, prerec_curve,general_scores = self._calculate_TPR_FPR_AUC(prediction_table)
            performance = {'AUC':auc_curve,'PreRec':prerec_curve,'general_scores':general_scores}
            self.Viraltools_benchmark_bin[prediction_tool+'_'+summary_metric] = performance
        else:
            for prediction_table in prediction_data_frames:
                prediction_tool = prediction_table[0][0]
                auc_curve, prerec_curve,general_scores = self._calculate_TPR_FPR_AUC(prediction_table)
                performance = {'AUC':auc_curve,'PreRec':prerec_curve,'general_scores':general_scores}

                # Save results
                self.Viraltools_benchmark_bin[prediction_tool+'_'+summary_metric] = performance


    @classmethod
    def _return_viraltool_prediction_tables(cls,tool,contigs,contig_type_lookup):
        
        viral_annotation_tools_cutoffs = {'virfinder':0.9,'virusseeker':0.9,'virsorter2':0.9,'deepvirfinder':0.9,'viralverify':15}

        predictions_table = []
        predicted_contigs = {'viral':set(),'bacterial':set()}
        for contig in contigs.values():
            contig_annotation = contig.contig_annotation
            if not tool in contig_annotation:           # Predicter did not recognise the contig as anything
                continue 
            viral_score = contig_annotation[tool][1]
            true_contig_label = contig_type_lookup[contig.name]
            if viral_score >= viral_annotation_tools_cutoffs[tool]:
                predicted_contigs['viral'].add(contig)
                predictions_table.append([tool,contig.name,true_contig_label,'viral',viral_score])
            else:
                predicted_contigs['bacterial'].add(contig)
                predictions_table.append([tool,contig.name,true_contig_label,'bacterial',viral_score])
        return predictions_table

    @classmethod
    def _confusion_matrix(cls,contig_types, predicted_contigs):
        '''ON CONTIG LEVEL: Given viral-predicted contigs for a tool, returns TP, TN, FP, FN
        '''
        true_viral_intersection = contig_types['viral'].intersection(predicted_contigs['viral'])
        true_bacterial_intersection = contig_types['bacterial'].intersection(predicted_contigs['bacterial'])
        actual_viral_breadth = sum([contig.len for contig in contig_types['viral']])
        actual_bacterial_breadth =  sum([contig.len for contig in contig_types['bacterial']])
        
        true_positives = sum([contig.len for contig in true_viral_intersection])
        false_positives = actual_viral_breadth - true_positives
        true_negatives =  sum([contig.len for contig in true_bacterial_intersection])
        false_negatives = actual_bacterial_breadth-true_negatives

        return true_positives, false_positives, true_negatives, false_negatives

    @classmethod
    def _calculate_TPR_FPR_AUC(cls,prediction_table):
        '''ON CONTIG LEVEL: Calculate AUC and Recall+Precision'''
        from sklearn.metrics import precision_recall_curve, roc_curve, auc, matthews_corrcoef

        prediction_tool = prediction_table[0][0]
        truth = [item[2] for item in prediction_table]
        predictions = [item[3] for item in prediction_table]
        scores = [item[4] for item in prediction_table]

        ### Transform Prediction labels
        binary_labels = _np.where(_np.array(truth) == 'bacterial',0,1)
        binary_predictions = _np.where(_np.array(predictions) == 'bacterial',0,1)

        fpr, tpr, roc_thresholds = roc_curve(binary_labels, scores)
        precision, recall, prerec_thresholds = precision_recall_curve(binary_labels,scores)
        precision, recall = precision[:-1],recall[:-1]
        auc_score = auc(fpr, tpr)
        mcc_score = matthews_corrcoef(truth,predictions)
        F1_score = f1_score(binary_labels,binary_predictions)
        AUCScores = _namedtuple('AUCScores', 'tool fpr tpr threshold')
        PreRecScores = _namedtuple('PreRecScores', 'tool precision recall threshold')
        general_scores = {'Tool':prediction_tool,'AUC':auc_score,'MCC':mcc_score,'F1':F1_score}

        auc_curve = []
        for i in range(len(fpr)):
            auc_curve.append(AUCScores(prediction_tool,round(fpr[i],3),round(tpr[i],3), roc_thresholds[i] ))

        prerec_curve = []
        for i in range(len(precision)):
            prerec_curve.append(PreRecScores(prediction_tool,round(precision[i],3), round(recall[i],3),prerec_thresholds[i]))

        return auc_curve, prerec_curve,general_scores
    
    @classmethod
    def _return_RF_prediction_dataframe(cls,genomes):
        '''ON GENOME LEVEL: Return sparse-matrix for Prediction with PHAMB model'''
        annotation_types = ['micompletehmm','voghmm','deepvirfinder']
        df = []
        genome_order = []
        for genome in genomes.values():
            genome_order += [genome.name]
            row = [genome.breadth]
            for type in annotation_types:
                value = 0
                if type in genome.genome_annotation:
                    if type == 'deepvirfinder':
                        value = genome.genome_annotation[type]['weighted_mean_score']
                    else:
                        value = genome.genome_annotation[type]
                row += [value]
            df.append(row)
        sparse_df = sparse.csr_matrix(df)
        return genome_order, sparse_df

    @classmethod
    def _run_RF_model(cls,RF_model,genometype_lookup,genome_order, sparse_df):
        '''ON GENOME LEVEL: Runs RF predictive-PHAMB model'''

        print('Loading Model and annotation table')
        trained_model = joblib.load(RF_model)

        predicted_genome_labels = trained_model.predict(sparse_df)
        prediction_probabilities = trained_model.predict_proba(sparse_df)
        predicted_genome_labels = [label.lower() for label in list(predicted_genome_labels) ] 
        rows = []
        for i, genome_name in enumerate(genome_order):
            rows.append(['PHAMB',genome_name, genometype_lookup[genome_name] ,predicted_genome_labels[i], prediction_probabilities[i][1]])
        return rows
    
    
    ### Code to return contig table to assert performance 
    @classmethod
    def _return_viraltool_prediction_tables_binlevel(cls, genomes, summary_metric):
        '''ON GENOME LEVEL: Returns a list of dataframes with Predicted Genome Labels for various Viral annotation tools 
        '''
        genometype_lookup = dict() 
        for genome in genomes.values():
            genometype_lookup[genome.name] = 'bacterial' if 'GCF' in genome.name else 'viral'

        viral_annotation_tools = ['virusseeker','virsorter2','virfinder','deepvirfinder','viralverify']
        viral_annotation_tools_cutoffs = {'virfinder':0.9,'virusseeker':0.9,'virsorter2':0.9,'deepvirfinder':0.9,'viralverify':15}
        data_frame_list = []
        for annotation_type in viral_annotation_tools:
            rows = []
            for genome in genomes.values():
                row = [annotation_type,genome.name, genometype_lookup[genome.name]  ]

                if not annotation_type in genome.genome_annotation:
                    row += ['bacterial',0,summary_metric]
                    rows.append(row)
                else:

                    ViralScores = genome.genome_annotation[annotation_type]     # Dictionary
                    bin_score = ViralScores[summary_metric]
                    if bin_score > viral_annotation_tools_cutoffs[annotation_type]:
                        row += ['viral',bin_score,summary_metric]
                    else:
                        row += ['bacterial',bin_score,summary_metric]
                    rows.append(row)
            data_frame_list.append(rows)
        
        return data_frame_list
    

    @staticmethod
    def _write_Viraltools_benchmark_results(results_directory,performance,bin_level=False):
        
        if not _os.path.exists(results_directory):
            _os.makedirs(results_directory)

        ### General scores
        file = _os.path.join(results_directory,'general_scores.tsv')

        with open(file,'w') as filehandle:
            filehandle.write('tool\tmetric\tscore\tsummary\n')
            for tool in performance:
                scores = performance[tool]['general_scores']
                if bin_level:
                    summary_metric = '_'.join(tool.split('_')[1:])
                    tool = tool.split('_')[0]
                else:
                    summary_metric = 'NA'
                for metric in scores:
                    if metric != 'Tool':
                        lineout = [tool,metric,round(scores[metric],2),summary_metric]
                        lineout = [str(i) for i in lineout]
                        filehandle.write('\t'.join(lineout)+'\n')
        
        ### AUC 
        file = _os.path.join(results_directory,'AUC.tsv')
        with open(file,'w') as filehandle:
            filehandle.write('tool\tfpr\ttpr\tthreshold\tsummary\n')
            for tool in performance:
                scores = performance[tool]['AUC']
                if bin_level:
                    summary_metric = '_'.join(tool.split('_')[1:])
                    tool = tool.split('_')[0]
                else:
                    summary_metric = 'NA'
                for tup in scores:
                    lineout = [tool,tup.fpr, tup.tpr, tup.threshold,summary_metric]
                    lineout = [str(i) for i in lineout]
                    filehandle.write('\t'.join(lineout)+'\n')

        ### Recall / Precision 
        file = _os.path.join(results_directory,'ROC.tsv')
        with open(file,'w') as filehandle:
            filehandle.write('tool\tprecision\trecall\tthreshold\tsummmary\n')
            for tool in performance:
                scores = performance[tool]['PreRec']
                if bin_level:
                    summary_metric = '_'.join(tool.split('_')[1:])
                    tool = tool.split('_')[0]
                else:
                    summary_metric = 'NA'
                for tup in scores:
                    lineout = [tool,tup.precision, tup.recall, tup.threshold,summary_metric]
                    lineout = [str(i) for i in lineout]
                    filehandle.write('\t'.join(lineout)+'\n')






class viral_annotation_tool_parsers:
    '''Set of methods for parssing the various output formats of viral predictors'''
    def __init__(self):
        pass
    
    @staticmethod
    def _parse_hmm_row(line):
        '''Both VOG and MiComplete HMM-searches - These are done on protein level'''
        items = line[:-1].split()
        protein_name, target, evalue, score= items[0], items[2],items[4], items[5]
        contig_name = protein_name.split('_')[0]
        score = round(float(score),2)
        evalue = float(evalue)
        annotation = (contig_name, target, score, evalue)
        if score >= 50:
            return annotation

    @staticmethod
    def _parse_virusseker_row(line):
        '''Will ONLY work with Virusseeker format: 3 columns with Contig-name\tprediction\tscore'''
        contig_name, prediction, score = line[:-1].split()
        score =round(float(score),2)
        if prediction != 'Phage':
            score = 0
        annotation = (contig_name, score, prediction)
        return annotation

    @staticmethod
    def _parse_vs2_row(line):
        '''Will ONLY work with Virsorter2 format'''
        contig_name, dsDNAphage,NCLDV, RNA, ssDNA, lavidaviridae, max_score, max_score_group, length, hallmark, viral, celluar = line[:-1].split()
        nhallmarks = int(hallmark)
        score = 0 if max_score == 'NaN' else max_score
        score =round(float(score),2)
        contig_name = contig_name.split('||')[0]
        annotation = (contig_name, score, nhallmarks)
        return annotation
    
    @staticmethod
    def _parse_virfinder_row(line):
        '''Will ONLY work with Virfinder format'''
        contig_name, length, score, pvalue = line[:-1].split()
        score =round(float(score),2)
        pvalue = float(pvalue)
        annotation = (contig_name, score, pvalue)
        return annotation
    
    @staticmethod
    def _parse_viralverify_row(line):
        '''Will ONLY work with Viralverify result format'''
        items = line.strip().split(',')
        if len(items) == 6:
            contig_name, prediction, length, circular, score, pfam_hits = items
        else:
             contig_name, prediction, length, score, pfam_hits = items
        if score == '-':
            score = 0
        elif score == '+':
            score = 0
        score =round(float(score),2)
        annotation = (contig_name, score, prediction)
        return annotation
        
    @staticmethod
    def _parse_dvf_row(line):
        '''Will ONLY work with Deepvirfinder result format'''
        contig_name, length, score, pvalue = line[:-1].split()
        score =round(float(score),2)
        pvalue = float(pvalue)
        annotation = (contig_name, score, pvalue)
        return annotation
    
    @staticmethod
    def _parse_vibrant_row(line):
        '''Will ONLY work with Vibrant result format'''
        contig_name, virustype, quality = line[:-1].split()
        annotation = (contig_name, virustype, quality)
        return annotation

    @staticmethod
    def _parse_checkv(line):
        '''Will ONLY work with CheckV output from v. X'''
        contig_name, contig_len, proviral_len, aai_expected, aai_completeness, aai_confidence, aai_error, aai_num_hits, aai_top_hit, aai_id, aai_af, hmm_upper, hmm_lower, hmm_num_hits, kmer_freq  = line[:-1].split()
        annotation = (contig_name, float(aai_completeness), float(aai_id), float(aai_af),aai_top_hit)
        return annotation


def checkv_genome_consensus(viral_annotation,fileout):
    '''Calculate the most common CheckV contig hit for a Genome'''
    genome_checkv_hits = dict()
    with open(fileout,'w') as out:
        out.write('Genome\ttop_hit\tproportion\tncontigs\n')
        for genome_name in viral_annotation.genomes:
            lineout = [genome_name]
            contig_hits = []
            for contig in viral_annotation.genomes[genome_name].contigs:
                top_hit = contig.contig_annotation['checkv'][-2]
                contig_hits += [top_hit]
            if len(contig_hits) == 0:
                genome_checkv_hits[genome_name] = ('NA','NA','NA')
                lineout += ['NA','NA','NA']
                out.write('\t'.join(lineout)+'\n')
            else:
                top_hit, ncontigs = _collections.Counter(contig_hits).most_common()[0]
                top_hit_proportion = ncontigs/len(contig_hits)
                genome_checkv_hits[genome_name] = (top_hit,top_hit_proportion,ncontigs)
                lineout += [top_hit,top_hit_proportion,ncontigs]
                out.write('\t'.join([str(i) for i in lineout])+'\n')
    return genome_checkv_hits
        

def filter_clusters(clusters, reference, minsize, mincontigs, checkpresence=True):
    """Creates a shallow copy of clusters, but without any clusters with a total size
    smaller than minsize, or fewer contigs than mincontigs.
    If checkpresence is True, raise error if a contig is not present in reference, else
    ignores it when counting cluster size.
    """

    filtered = dict()
    for binname, contignames in clusters.items():
        if len(contignames) < mincontigs:
            continue
        size = 0
        for contigname in contignames:
            contig = reference.contigs.get(contigname)
            if contig is not None:
                size += len(contig)
            elif checkpresence:
                raise KeyError('Contigname {} not in reference'.format(contigname))
            else:
                pass
        if size >= minsize:
            filtered[binname] = contignames.copy()

    return filtered

def write_concat_bins(directory, bins, fastadict, compressed=False, maxbins=250, minsize=5000): 
    """Writes bins as FASTA files in a directory, one file per bin.
    Inputs:
        directory: Directory to create or put files in
        bins: {'name': {set of contignames}} dictionary (can be loaded from
        clusters.tsv using vamb.cluster.read_clusters)
        fastadict: {contigname: FastaEntry} dict as made by `loadfasta`
        compressed: Sequences in dict are compressed [False]
        maxbins: None or else raise an error if trying to make more bins than this [250]
        minsize: Minimum number of nucleotides in cluster to be output [0]
    Output: None
    """
    import os as _os
    import gzip as _gzip
    from _vambtools import FastaEntry
    import random

    # Safety measure so someone doesn't accidentally make 50000 tiny bins
    # If you do this on a compute cluster it can grind the entire cluster to
    # a halt and piss people off like you wouldn't believe.
    if maxbins is not None and len(bins) > maxbins:
        raise ValueError('{} bins exceed maxbins of {}'.format(len(bins), maxbins))

    # Check that the directory is not a non-directory file,
    # and that its parent directory indeed exists
    abspath = _os.path.abspath(directory)
    parentdir = _os.path.dirname(abspath)

    if parentdir != '' and not _os.path.isdir(parentdir):
        raise NotADirectoryError(parentdir)

    if _os.path.isfile(abspath):
        raise NotADirectoryError(abspath)

    if minsize < 0:
        raise ValueError("Minsize must be nonnegative")

    # Check that all contigs in all bins are in the fastadict
    allcontigs = set()

    for contigs in bins.values():
        allcontigs.update(set(contigs))

    allcontigs -= fastadict.keys()
    if allcontigs:
        nmissing = len(allcontigs)
        raise IndexError('{} contigs in bins missing from fastadict'.format(nmissing))

    # Make the directory if it does not exist - if it does, do nothing
    try:
        _os.mkdir(directory)
    except FileExistsError:
        pass
    except:
        raise
    
    bins_entries = []
    # Now actually print all the contigs to files
    for binname, contigs in bins.items():
        
        # Concatenate sequences of the bin
        concat_sequence = bytearray()
        for contig in contigs:
            entry = fastadict[contig]
            if compressed:
                uncompressed = bytearray(_gzip.decompress(entry.sequence))
                concat_sequence += uncompressed
            else:
                uncompressed = bytearray(entry.sequence)
                concat_sequence += uncompressed

        bin_entry = FastaEntry(binname, concat_sequence)           
        # Skip bin if it's too small
        if len(bin_entry.sequence) < minsize:
            continue
        bins_entries.append(bin_entry)
    
    random.shuffle(bins_entries)
    print('Writing:',len(bins_entries) ,'bins to file')
    filename = _os.path.join(directory, 'vamb_bins.1.fna')
    i = 1
    j = 1
    file = open(filename,'w')
    for entry in bins_entries:
        if i % 100000 == 0:
            j += 1 
            file.close()
            filename = _os.path.join(directory, 'vamb_bins.' + str(j) + '.fna')
            file = open(filename,'w')
        i += 1
        print(entry.format(),file=file)
