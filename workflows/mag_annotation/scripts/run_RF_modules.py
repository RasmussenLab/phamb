#!/bin/python
'''Helper modules'''
import vambtools as _vambtools
import collections as _collections
import os
import numpy as _np

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
    __slots__ = ['name', 'breadth', 'contigs','totalsize','genome_annotation']

    def __init__(self, name):
        self.name = name
        self.contigs = set()
        self.breadth = 0
        self.totalsize = 0
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
    
    def gettotalsize(self):
        "This calculates the total size of Genome (potentially redundant)"
        bysubject = _collections.defaultdict(list)
        for contig in self.contigs:
            bysubject[contig.subject].append(contig)
        genomesize = 0
        for contiglist in bysubject.values():
            for contig in contiglist:
                genomesize += contig.end
        self.totalsize = genomesize

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
        for genome in genomes:
            genome.gettotalsize()

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


class Viral_annotation:
    '''A module for parsing all kinds of annotation pr. contig for a Reference Genome or VAMB bin'''
    # Instantiate with any iterable of Genomes
    def __init__(self, genomes, annotation_files):
        '''genomes is either a Reference Genome or VAMB bin'''

        self.contigs = genomes.contigs      # Make shallow copy of contigs 
        self.genomes = genomes.genomes      # Make shallow copy of genomes

        ### Check files
        for f in annotation_files.values():
            if not os.path.exists(f):                
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
                    if annotation_tuple[0] in self.contigs:
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


class viral_annotation_tool_parsers:
    '''Set of methods for parssing the various output formats of viral predictors'''
    def __init__(self):
        pass
    
    @staticmethod
    def _parse_hmm_row(line):
        '''Both VOG and MiComplete HMM-searches and --tblout files '''
        items = line[:-1].split()
        protein_name, target, evalue, score= items[0], items[2],items[4], items[5]
        contig_name = '_'.join(protein_name.split('_')[:-1]) # Remove protein number 
        score = round(float(score),2)
        evalue = float(evalue)
        annotation = (contig_name, target, score, evalue)
        if score >= 50:
            return annotation

    @staticmethod
    def _parse_virusseker_row(line):
        '''Will ONLY work with Virusseeker format: 3 columns with '''
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
        contig_name, length, score, pvalue = line[:-1].split("\t")
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


