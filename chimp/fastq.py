from Bio import SeqIO
from collections import defaultdict
import gzip
from model.fastq import FastqRead, FastqFiles
from model.tile import Tile
import os
import progressbar
import pysam
import subprocess


def parse_fastq_gzip_file(path):
    with gzip.open(path) as f:
        for dna in SeqIO.parse(f, 'fastq'):
            yield FastqRead(dna)
    del f


def load_tiles(fastq_files):
    """
    We need to do this, but not using the raw data. We should use the parsed files so we can get sequence info too

    """
    reads = defaultdict(list)
    for filename in fastq_files:
        for r in parse_fastq_gzip_file(filename):
            reads[r.region].append(r)
    for region, fastq_reads in reads.items():
        yield Tile(region, fastq_reads)


class FastqReadClassifier(object):
    def __init__(self, bowtie_path):
        clean_path = bowtie_path.rstrip(os.path.sep)
        self.name = os.path.basename(clean_path)
        self._common_command = ('bowtie2', '--local', '-p 15', '--no-unal', '-x %s' % clean_path)

    def paired_call(self, fastq_file_1, fastq_file_2):
        command = self._common_command + ('-1 ' + fastq_file_1,
                                          '-2 ' + fastq_file_2,
                                          '-S /tmp/chimp.sam',
                                          '2>&1 | tee /tmp/error.txt')
        return self._run(command)

    def single_call(self, fastq_file):
        command = self._common_command + ('-U ' + fastq_file,)
        return self._run(command)

    def _run(self, command):
        subprocess.call(' '.join(command), shell=True)
        sam_command = 'samtools view -bS /tmp/chimp.sam | samtools sort - /tmp/final'
        subprocess.call(sam_command, shell=True)
        subprocess.call('samtools index /tmp/final.bam', shell=True)
        for r in pysam.Samfile('/tmp/final.bam'):
            yield r.qname


def classify_fastq_reads(classifier_path, fastq_files):
    classifier = FastqReadClassifier(classifier_path)
    reads = set()
    pbar = progressbar.ProgressBar(max_value=fastq_files.alignment_length)
    current = 1
    pbar.update(current)
    for file1, file2 in fastq_files.paired:
        for read in classifier.paired_call(file1, file2):
            reads.add(read)
        current += 1
        pbar.update(current)
    for file1 in fastq_files.single:
        for read in classifier.single_call(file1):
            reads.add(read)
        current += 1
        pbar.update(current)
    pbar.update(current)
    return classifier.name, reads


def classify_all_reads(classifier_paths, fastq_files):
    classified_reads = {}
    for path in classifier_paths:
        name, reads = classify_fastq_reads(path, fastq_files)
        assert name != 'unclassified', '"unclassified" cannot be used as a fastq read classifier name'
        classified_reads[name] = reads
    return classified_reads


def stream_all_read_names(fastq_files):
    # streams sets of read names for each fastq file
    for filename in fastq_files:
        yield {r.name for r in parse_fastq_gzip_file(filename)}


def save_classified_reads(name, reads):
    with open('classified_reads/%s' % name, 'w+') as f:
        for read in reads:
            f.write('%s\n' % read)


def main(fastq_directory, classifier_paths):
    """
    Parses fastq files and creates text files containing read names that belong to each source
    Currently this is just used to find phiX

    """
    filenames = [os.path.join(fastq_directory, filename) for filename in os.listdir(fastq_directory)]
    fastq_files = FastqFiles(filenames)
    classified_reads = classify_all_reads(classifier_paths, fastq_files)

    all_classified_reads = set()
    for name, reads in classified_reads.items():
        # write the reads to disk for later use
        save_classified_reads(name, reads)
        # keep a single set with all the reads so that we can figure out which reads aren't classified
        all_classified_reads.update(reads)

    # now figure out which reads remain unclassified and save them to a catch-all bucket
    all_unclassified_reads = set()
    for all_read_names in stream_all_read_names(fastq_files):
        all_unclassified_reads.update(all_read_names.difference(all_classified_reads))
    save_classified_reads('unclassified', all_unclassified_reads)


if __name__ == '__main__':
    main('SA15243/all_fastqs', ['/home/jim/Desktop/ngs/phix/bowtie2/phix'])
