import gzip
from Bio import SeqIO
from model.fastq import FastqRead


def parse_fastq_gzip_file(path):
    with gzip.open(path) as f:
        for dna in SeqIO.parse(f, 'fastq'):
            yield FastqRead(dna.name, str(dna.seq))
