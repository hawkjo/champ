from Bio import SeqIO
from collections import defaultdict
import gzip
from model.fastq import FastqRead, FastqFiles
from model.tile import Tile
import os


def parse_fastq_lines(f):
    for dna in SeqIO.parse(f, 'fastq'):
        yield FastqRead(dna)


def parse_fastq_gzip_file(path):
    with gzip.open(path) as f:
        for fastq_read in parse_fastq_lines(f):
            yield fastq_read
    del f


def load_fastq_filenames(directory):
    filenames = [os.path.join(directory, filename) for filename in os.listdir(directory)]
    sizes = map(os.path.getsize, filenames)
    data = {filename: size for filename, size in zip(filenames, sizes)}
    return FastqFiles(data)


def load_all_fastq_data(directory):
    filenames = load_fastq_filenames(directory)
    reads = defaultdict(list)
    for n, filename in enumerate(filenames):
        for r in parse_fastq_gzip_file(filename):
            reads[r.region].append(r)
    return reads


def load_tiles(directory):
    reads = load_all_fastq_data(directory)
    for region, fastq_reads in reads.items():
        yield Tile(region, fastq_reads)
