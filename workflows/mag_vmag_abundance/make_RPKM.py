# This script calculates RPKM, when paired end reads are mapped to a contig
# catalogue with BWA MEM. It will not be accurate with single end reads or
# any other mapper than BWA MEM.

# Theory:
# We want a simple way to estimate abundance of redundant contig catalogues.
# Earlier we used to run analysis on deduplicated gene catalogues, but since
# both depth and kmer composition are only stable for longer contigs, we have
# moved to contig catalogues. We have not found a way of deduplicating contigs.

# For this we have until now used two methods:
# 1) Only counting the primary hits. In this case the read will never be
# assigned to any contig which differ by just 1 basepair. Even for
# identical contigs, reads are assigned randomly which causes noise.

# 2) Using MetaBAT's jgi_summarize_bam_contig_depths, a script which is not
# documented and we cannot figure out how works. When testing with small
# toy data, it produces absurd results.

# This script is an attempt to take an approach as simple as possible while
# still being sound technically. We simply count the number of reads in a
# contig normalized by contig length and total number of reads.

# We look at all hits, including secondary hits. We do not discount partial
# alignments. Also, if a read maps to N contigs, we count each hit as 1/N reads.
# The reason for all these decisions is that if the aligner believes it's a hit,
# we believe the contig is present.

# We do not take varying insert sizes into account. It is unlikely that
# any contig with enough reads to provide a reliable estimate of depth would,
# by chance, only recruit read pairs with short or long insert size. So this
# will average out over all contigs.

# We count each read independently, because BWA MEM often assigns mating reads
# to different contigs..

__doc__ = """Estimate RPKM (depths) from BAM files of reads mapped to contigs.
Usage:
>>> bampaths = ['/path/to/bam1.bam', '/path/to/bam2.bam', '/path/to/bam3.bam']
>>> rpkms = read_bamfiles(bampaths)
"""

import pysam as _pysam
import shutil
import sys as _sys
import os as _os
import multiprocessing as _multiprocessing
import numpy as _np
import time as _time
from hashlib import md5 as _md5
import argparse
import vamb.vambtools as _vambtools
import vamb


parser = argparse.ArgumentParser(description='''
    Annotating contigs to filter off possible non viral sequences
   ''')
parser.add_argument('-b', help='List of BAM files')
parser.add_argument('-f', help='Reference fasta')
parser.add_argument('-o', help='Directory out')
parser.add_argument('-t', help='Threads')


DEFAULT_SUBPROCESSES = 28


def read_contigs(filehandle, minlength=2000):
    """Parses a FASTA file open in binary reading mode.
    Input:
        filehandle: Filehandle open in binary mode of a FASTA file
        minlength: Ignore any references shorter than N bases [100]
    Outputs:
        tnfs: An (n_FASTA_entries x 103) matrix of tetranucleotide freq.
        contignames: A list of contig headers
        lengths: A Numpy array of contig lengths
    """

    if minlength < 4:
        raise ValueError('Minlength must be at least 4, not {}'.format(minlength))

    raw = _vambtools.PushArray(_np.float32)
    projected = _vambtools.PushArray(_np.float32)
    lengths = _vambtools.PushArray(_np.int)
    contignames = list()

    entries = _vambtools.byte_iterfasta(filehandle)

    for entry in entries:
        if len(entry) < minlength:
            continue

        lengths.append(len(entry))
        contignames.append(entry.header)

    # Don't use reshape since it creates a new array object with shared memory
    lengths_arr = lengths.take()

    return contignames, lengths_arr

def mergecolumns(pathlist):
    """Merges multiple npz files with columns to a matrix.
    All paths must be npz arrays with the array saved as name 'arr_0',
    and with the same length.
    Input: pathlist: List of paths to find .npz files to merge
    Output: Matrix with one column per npz file
    """

    if len(pathlist) == 0:
        return _np.array([], dtype=_np.float32)

    for path in pathlist:
        if not _os.path.exists(path):
            raise FileNotFoundError(path)

    first = _np.load(pathlist[0])['arr_0']
    length = len(first)
    ncolumns = len(pathlist)

    result = _np.zeros((length, ncolumns), dtype=_np.float32)
    result[:,0] = first

    for columnno, path in enumerate(pathlist[1:]):
        column = _np.load(path)['arr_0']
        if len(column) != length:
            raise ValueError("Length of data at {} is not equal to that of "
                             "{}".format(path, pathlist[0]))
        result[:,columnno + 1] = column

    return result

def _identity(segment):
    "Return the nucleotide identity of the given aligned segment."
    mismatches, matches = 0, 0
    for kind, number in segment.cigartuples:
        # 0, 7, 8, is match/mismatch, match, mismatch, respectively
        if kind in (0, 7, 8):
            matches += number
        # 1, 2 is insersion, deletion
        elif kind in (1, 2):
            mismatches += number
    matches -= segment.get_tag('NM')
    return matches / (matches+mismatches)

def _filter_segments(segmentiterator, minscore, minid):
    """Returns an iterator of AlignedSegment filtered for reads with low
    alignment score.
    """

    for alignedsegment in segmentiterator:
        # Skip if unaligned or suppl. aligment
        if alignedsegment.flag & 0x804 != 0:
            continue

        if minscore is not None and alignedsegment.get_tag('AS') < minscore:
            continue

        if minid is not None and _identity(alignedsegment) < minid:
            continue

        yield alignedsegment

def count_reads(bamfile, minscore=None, minid=None):
    """Count number of reads mapping to each reference in a bamfile,
    optionally filtering for score and minimum id.
    Multi-mapping reads MUST be consecutive in file, and their counts are
    split among the references.
    Inputs:
        bamfile: Open pysam.AlignmentFile
        minscore: Minimum alignment score (AS field) to consider [None]
        minid: Discard any reads with ID lower than this [None]
    Output: Float32 Numpy array of read counts for each reference in file.
    """
    # Use 64-bit floats for better precision when counting
    readcounts = _np.zeros(len(bamfile.lengths))

    # Initialize with first aligned read - return immediately if the file
    # is empty
    filtered_segments = _filter_segments(bamfile, minscore, minid)
    try:
        segment = next(filtered_segments)
        read_name = segment.query_name
        multimap = 1.0
        reference_ids = [segment.reference_id]
    except StopIteration:
        return readcounts.astype(_np.float32)

    # Now count up each read in the BAM file
    for segment in filtered_segments:
        # If we reach a new read_name, we tally up the previous read
        # towards all its references, split evenly.
        if segment.query_name != read_name:
            read_name = segment.query_name
            to_add = 1.0 / multimap
            for reference_id in reference_ids:
                readcounts[reference_id] += to_add
            reference_ids.clear()
            multimap = 0.0

        multimap += 1.0
        reference_ids.append(segment.reference_id)

    # Add final read
    to_add = 1.0 / multimap
    for reference_id in reference_ids:
        readcounts[reference_id] += to_add

    return readcounts.astype(_np.float32)

def calc_rpkm(counts, lengths, minlength=None):
    """Calculate RPKM based on read counts and sequence lengths.
    Inputs:
        counts: Numpy vector of read counts from count_reads
        lengths: Iterable of contig lengths in same order as counts
        minlength [None]: Discard any references shorter than N bases
    Output: Float32 Numpy vector of RPKM for all seqs with length >= minlength
    """
    lengtharray = _np.array(lengths)
    if len(counts) != len(lengtharray):
        raise ValueError("counts length and lengths length must be same")

    millionmappedreads = counts.sum() / 1e6

    # Prevent division by zero
    if millionmappedreads == 0:
        rpkm = _np.zeros(len(lengtharray), dtype=_np.float32)
    else:
        kilobases = lengtharray / 1000
        rpkm = (counts / (kilobases * millionmappedreads)).astype(_np.float32)

    # Now filter away small contigs
    if minlength is not None:
        lengthmask = lengtharray >= minlength
        rpkm = rpkm[lengthmask]

    return rpkm

def _hash_refnames(refnames):
    "Hashes an iterable of strings of reference names using MD5."
    hasher = _md5()
    for refname in refnames:
        hasher.update(refname.encode().rstrip())

    return hasher.digest()

def _check_bamfile(path, bamfile, refhash, minlength):
    "Checks bam file for correctness (refhash and sort order). To be used before parsing."
    # If refhash is set, check ref hash matches what is found.
    if refhash is not None:
        if minlength is None:
            refnames = bamfile.references
        else:
            pairs = zip(bamfile.references, bamfile.lengths)
            refnames = (ref for (ref, len) in pairs if len >= minlength)

        hash = _hash_refnames(refnames)
        if hash != refhash:
            errormsg = ('BAM file {} has reference hash {}, expected {}. '
                        'Verify that all BAM headers and FASTA headers are '
                        'identical and in the same order.')
            raise ValueError(errormsg.format(path, hash.hex(), refhash.hex()))

    # Check that file is unsorted or sorted by read name.
    hd_header = bamfile.header.get("HD", dict())
    sort_order = hd_header.get("SO")
    if sort_order in ("coordinate", "unknown"):
        errormsg = ("BAM file {} is marked with sort order '{}', must be "
                    "unsorted or sorted by readname.")
        raise ValueError(errormsg.format(path, sort_order))


def _get_contig_rpkms(inpath, outpath, refhash, minscore, minlength, minid):
    """Returns  RPKM (reads per kilobase per million mapped reads)
    for all contigs present in BAM header.
    Inputs:
        inpath: Path to BAM file
        outpath: Path to dump depths array to or None
        refhash: Expected reference hash (None = no check)
        minscore: Minimum alignment score (AS field) to consider
        minlength: Discard any references shorter than N bases
        minid: Discard any reads with ID lower than this
    Outputs:
        path: Same as input path
        rpkms:
            If outpath is not None: None
            Else: A float32-array with RPKM for each contig in BAM header
        length: Length of rpkms array
    """

    bamfile = _pysam.AlignmentFile(inpath, "rb")
    _check_bamfile(inpath, bamfile, refhash, minlength)
    counts = count_reads(bamfile, minscore, minid)
    rpkms = calc_rpkm(counts, bamfile.lengths, minlength)
    bamfile.close()

    # If dump to disk, array returned is None instead of rpkm array
    if outpath is not None:
        arrayresult = None
        _np.savez_compressed(outpath, rpkms)
    else:
        arrayresult = rpkms

    return inpath, arrayresult, len(rpkms)

def read_bamfiles(paths, dumpdirectory=None, refhash=None, minscore=None, minlength=None,
                  minid=None, subprocesses=DEFAULT_SUBPROCESSES, logfile=None):
    "Placeholder docstring - replaced after this func definition"
    
    paths_list = []
    with open(paths,'r') as infile:
        for line in infile:
            paths_list.append(line.strip())

    sample_order = _np.array([ _os.path.basename(i).replace('.bam','') for i in paths_list ])

    paths=paths_list
    # Define callback function depending on whether a logfile exists or not
    if logfile is not None:
        def _callback(result):
            path, rpkms, length = result
            print('\tProcessed', path, file=logfile)
            logfile.flush()

    else:
        def _callback(result):
            pass

    # Bam files must be unique.
    if len(paths) != len(set(paths)):
        raise ValueError('All paths to BAM files must be unique.')

    # Bam files must exist
    for path in paths:
        if not _os.path.isfile(path):
            raise FileNotFoundError(path)

        # Create directory to dump in
        if not _os.path.exists(dumpdirectory):
            _os.mkdir(dumpdirectory)

    # Spawn independent processes to calculate RPKM for each of the BAM files
    processresults = list()

    # Queue all the processes
    with _multiprocessing.Pool(processes=subprocesses) as pool:
        for pathnumber, path in enumerate(paths):
            if dumpdirectory is None:
                outpath = None
            else:
                outpath = _os.path.join(dumpdirectory, str(pathnumber) + '.npz')

            arguments = (path, outpath, refhash, minscore, minlength, minid)
            processresults.append(pool.apply_async(_get_contig_rpkms, arguments,
                                                   callback=_callback))

        all_done, any_fail = False, False
        while not (all_done or any_fail):
            _time.sleep(5)
            all_done = all(process.ready() and process.successful() for process in processresults)
            any_fail = any(process.ready() and not process.successful() for process in processresults)

            if all_done:
                pool.close() # exit gently
            if any_fail:
                pool.terminate() # exit less gently

        # Wait for all processes to be cleaned up
        pool.join()

    # Raise the error if one of them failed.
    for path, process in zip(paths, processresults):
        if process.ready() and not process.successful():
            print('\tERROR WHEN PROCESSING:', path, file=logfile)
            print('Vamb aborted due to error in subprocess. See stacktrace for source of exception.')
            if logfile is not None:
                logfile.flush()
            process.get()

    ncontigs = None
    for processresult in processresults:
        path, rpkm, length = processresult.get()

        # Verify length of contigs are same for all BAM files
        if ncontigs is None:
            ncontigs = length
        elif length != ncontigs:
            raise ValueError('First BAM file has {} headers, {} has {}.'.format(
                             ncontigs, path, length))

    # If we did not dump to disk, load directly from process results to
    # one big matrix...
    if dumpdirectory is None:
        columnof = {p:i for i, p in enumerate(paths)}
        rpkms = _np.zeros((ncontigs, len(paths)), dtype=_np.float32)

        for processresult in processresults:
            path, rpkm, length = processresult.get()
            rpkms[:, columnof[path]] = rpkm

    # If we did, instead merge them from the disk
    else:
        dumppaths = [_os.path.join(dumpdirectory, str(i) + '.npz') for i in range(len(paths))]
        rpkms = mergecolumns(dumppaths)

    return rpkms, sample_order



if __name__ == "__main__":
    args = parser.parse_args()
    DEFAULT_SUBPROCESSES = args.t 

    with vamb.vambtools.Reader(args.f, 'rb') as tnffile:
            contignames, lengths_arr = read_contigs(tnffile, minlength=2000)

    with open('vamb_log_97.txt', 'w') as logfile:
        refhash = vamb.parsebam._hash_refnames(contignames)

        dumpdirectory = 'vamb_tmp_id_97'
        rpkms, sample_order = read_bamfiles(args.b,dumpdirectory=dumpdirectory, refhash=refhash,   minscore=None,minlength=2000,minid=0.97,subprocesses=DEFAULT_SUBPROCESSES,logfile=logfile)
        
        ### Save Contignames, Lengths and RPKM matrix.
        vamb.vambtools.write_npz(_os.path.join(args.o, 'contigs.npz'), contignames)
        vamb.vambtools.write_npz(_os.path.join(args.o, 'lengths.npz'), lengths_arr)
        vamb.vambtools.write_npz(_os.path.join(args.o, 'samples.npz'), sample_order)
        vamb.vambtools.write_npz(_os.path.join(args.o, 'RPKM.npz'), rpkms)
        
        shutil.rmtree(dumpdirectory)



class Args:
    def __init__(self,b,v,o,f):
        self.b =  b
        self.v = v 
        self.o = o
        self.f = f

args = Args('mapping_files.txt','g','05_binning/vamb_on_jgi_v3/HMP2','combined_assemblies/HMP2.fna')
