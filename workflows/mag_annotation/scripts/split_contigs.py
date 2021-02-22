
#!/bin/python
from Bio import SeqIO
import os
import gzip
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    required = parser.add_argument_group('Required arguments')
    
    required.add_argument('-c',
                          '--contigs_file',
                          metavar='',
                          required=True,
                          type=str,
                          help='Fasta file with contigs')
    (args, extra_args) = parser.parse_known_args()

    return args

### Plenty of near and efficient classes and functions from VAMBtools


def byte_iterfasta(filehandle, comment=b'#'):
    """Yields FastaEntries from a binary opened fasta file.
    Usage:
    >>> with Reader('/dir/fasta.fna', 'rb') as filehandle:
    ...     entries = byte_iterfasta(filehandle) # a generator
    Inputs:
        filehandle: Any iterator of binary lines of a FASTA file
        comment: Ignore lines beginning with any whitespace + comment
    Output: Generator of FastaEntry-objects from file
    """

    # Make it work for persistent iterators, e.g. lists
    line_iterator = iter(filehandle)
    # Skip to first header
    try:
        for probeline in line_iterator:
            stripped = probeline.lstrip()
            if stripped.startswith(comment):
                pass

            elif probeline[0:1] == b'>':
                break

            else:
                raise ValueError('First non-comment line is not a Fasta header')

        else: # no break
            raise ValueError('Empty or outcommented file')

    except TypeError:
        errormsg = 'First line does not contain bytes. Are you reading file in binary mode?'
        raise TypeError(errormsg) from None

    header = probeline[1:-1].decode()
    buffer = list()

    # Iterate over lines
    for line in line_iterator:
        if line.startswith(comment):
            pass

        elif line.startswith(b'>'):
            yield FastaEntry(header, bytearray().join(buffer))
            buffer.clear()
            header = line[1:-1].decode()

        else:
            buffer.append(line)

    yield FastaEntry(header, bytearray().join(buffer))

def read_contigs(filehandle, minlength=2000):
    """Parses a FASTA file open in binary reading mode.
    Input:
        filehandle: Filehandle open in binary mode of a FASTA file
        minlength: Ignore any references shorter than N bases [100]
    Outputs:
        contignames: A list of contig headers
        lengths: A Numpy array of contig lengths
    """

    if minlength < 4:
        raise ValueError('Minlength must be at least 4, not {}'.format(minlength))

    lengths = PushArray(np.int)
    contignames = list()

    entries = byte_iterfasta(filehandle)

    for entry in entries:
        if len(entry) < minlength:
            continue

        lengths.append(len(entry))
        contignames.append(entry.header)

    # Don't use reshape since it creates a new array object with shared memory
    lengths_arr = lengths.take()

    return contignames, lengths_arr



def write_npz(file, array):
    """Writes a Numpy array to an open file or path in .npz format
    Inputs:
        file: Open file or path to file
        array: Numpy array
    Output: None
    """
    np.savez_compressed(file, array)

class Reader:
    """Use this instead of `open` to open files which are either plain text,
    gzipped, bzip2'd or zipped with LZMA.
    Usage:
    >>> with Reader(file, readmode) as file: # by default textmode
    >>>     print(next(file))
    TEST LINE
    """

    def __init__(self, filename, readmode='r'):
        if readmode not in ('r', 'rt', 'rb'):
            raise ValueError("the Reader cannot write, set mode to 'r' or 'rb'")
        if readmode == 'r':
            self.readmode = 'rt'
        else:
            self.readmode = readmode

        self.filename = filename

        with open(self.filename, 'rb') as f:
            signature = f.peek(8)[:8]

        # Gzipped files begin with the two bytes 0x1F8B
        if tuple(signature[:2]) == (0x1F, 0x8B):
            self.filehandle = gzip.open(self.filename, self.readmode)

        # Else we assume it's a text file.
        else:
            self.filehandle = open(self.filename, self.readmode)

    def close(self):
        self.filehandle.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __iter__(self):
        return self.filehandle

class FastaEntry:
    """One single FASTA entry. Instantiate with string header and bytearray
    sequence."""

    basemask = bytearray.maketrans(b'acgtuUswkmyrbdhvnSWKMYRBDHV',
                                   b'ACGTTTNNNNNNNNNNNNNNNNNNNNN')
    __slots__ = ['header', 'sequence']

    def __init__(self, header, sequence):
        if len(header) > 0 and (header[0] in ('>', '#') or header[0].isspace()):
            raise ValueError('Header cannot begin with #, > or whitespace')
        if '\t' in header:
            raise ValueError('Header cannot contain a tab')

        masked = sequence.translate(self.basemask, b' \t\n\r')
        stripped = masked.translate(None, b'ACGTN')
        if len(stripped) > 0:
            bad_character = chr(stripped[0])
            msg = "Non-IUPAC DNA byte in sequence {}: '{}'"
            raise ValueError(msg.format(header, bad_character))

        self.header = header
        self.sequence = masked

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return '>{}\n{}'.format(self.header, self.sequence.decode())

    def format(self, width=60):
        sixtymers = range(0, len(self.sequence), width)
        spacedseq = '\n'.join([self.sequence[i: i+width].decode() for i in sixtymers])
        return '>{}\n{}'.format(self.header, spacedseq)

    def __getitem__(self, index):
        return self.sequence[index]

    def __repr__(self):
        return '<FastaEntry {}>'.format(self.header)


class PushArray:
    """Data structure that allows efficient appending and extending a 1D Numpy array.
    Intended to strike a balance between not resizing too often (which is slow), and
    not allocating too much at a time (which is memory inefficient).
    Usage:
    >>> arr = PushArray(numpy.float64)
    >>> arr.append(5.0)
    >>> arr.extend(numpy.linspace(4, 3, 3))
    >>> arr.take() # return underlying Numpy array
    array([5. , 4. , 3.5, 3. ])
    """

    __slots__ = ['data', 'capacity', 'length']

    def __init__(self, dtype, start_capacity=1<<16):
        self.capacity = start_capacity
        self.data = np.empty(self.capacity, dtype=dtype)
        self.length = 0

    def __len__(self):
        return self.length

    def _setcapacity(self, n):
        self.data.resize(n, refcheck=False)
        self.capacity = n

    def _grow(self, mingrowth):
        """Grow capacity by power of two between 1/8 and 1/4 of current capacity, though at
        least mingrowth"""
        growth = max(int(self.capacity * 0.125), mingrowth)
        nextpow2 = 1 << (growth - 1).bit_length()
        self._setcapacity(self.capacity + nextpow2)

    def append(self, value):
        if self.length == self.capacity:
            self._grow(64)

        self.data[self.length] = value
        self.length += 1

    def extend(self, values):
        lenv = len(values)
        if self.length + lenv > self.capacity:
            self._grow(lenv)

        self.data[self.length:self.length+lenv] = values
        self.length += lenv

    def take(self):
        "Return the underlying array"
        self._setcapacity(self.length)
        return self.data

    def clear(self, force=False):
        "Empties the PushArray. If force is true, also truncates the underlying memory."
        self.length = 0
        if force:
            self._setcapacity(0)
    



### Generic function to split Contigs file into sample-seperated files with contigs
def splitcontigs_to_samples(fasta_file,directory_out):
    
    current_sample = None
    outhandle = None
    samples = set()
    for record in SeqIO.parse(gzip.open(fasta_file, 'rt'), 'fasta'):
            record.description = record.id
            contig = record.id
            sample = contig.split('C')[0]
            samples.add(sample)

            if sample != current_sample:

                if not outhandle is None:
                    outhandle.close() 
                os.makedirs(os.path.join(directory_out,sample),exist_ok=True)
                outhandle = open('{}/{}/{}.fna'.format(directory_out,sample,sample),  'w')
                SeqIO.write(record, outhandle, 'fasta') 
                current_sample = sample
            else:
                SeqIO.write(record, outhandle, 'fasta') 
                current_sample = sample
    
    with open('sample_table.txt','w') as out:
        for s in samples:
            out.write(s+'\n')


if __name__ == "__main__":
    args = parse_arguments()
    splitcontigs_to_samples('contigs.fna.gz','assembly')



    with Reader('contigs.fna.gz', 'rb') as fastafile:
            contignames, lengths_arr = read_contigs(fastafile, minlength=2000)

    write_npz('contigs.npz',contignames)
    write_npz('contig_lengths.npz',lengths_arr)
