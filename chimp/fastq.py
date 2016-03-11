from Bio import SeqIO
from chimp.model.fastq import FastqRead
from collections import defaultdict, namedtuple
import gzip
import logging
import os
from progressbar import ProgressBar, Percentage, Counter, Bar
import pysam
import subprocess


log = logging.getLogger(__name__)


FastqName = namedtuple('FastqRead', ['name'])


def parse_fastq_gzip_file(path):
    with gzip.open(path) as f:
        for dna in SeqIO.parse(f, 'fastq'):
            yield FastqRead(dna)
    del f


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
        with open('/dev/null', 'w+') as devnull:
            kwargs = dict(shell=True, stderr=devnull, stdout=devnull)
            subprocess.call(' '.join(command), **kwargs)
            sam_command = 'samtools view -bS /tmp/chimp.sam | samtools sort - /tmp/final'
            subprocess.call(sam_command, **kwargs)
            subprocess.call('samtools index /tmp/final.bam', **kwargs)
            for r in pysam.Samfile('/tmp/final.bam'):
                yield r.qname


def classify_fastq_reads(classifier_path, fastq_files):
    progress_bar = get_progress_bar(fastq_files.alignment_length)
    classifier = FastqReadClassifier(classifier_path)
    log.info('Searching for reads for %s. This will take a while!' % classifier.name)

    current = 0
    reads = set()
    for file1, file2 in fastq_files.paired:
        for read in classifier.paired_call(file1, file2):
            reads.add(read)
        current += 1
        progress_bar.update(current)
    for file1 in fastq_files.single:
        for read in classifier.single_call(file1):
            reads.add(read)
        current += 1
        progress_bar.update(current)
    progress_bar.finish()
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


def load_unclassified_reads(fastq_files, all_classified_reads):
    all_unclassified_reads = set()
    progress_bar = get_progress_bar(len(fastq_files))
    for all_read_names in progress_bar(stream_all_read_names(fastq_files)):
        all_unclassified_reads.update(all_read_names.difference(all_classified_reads))
    return all_unclassified_reads


def get_progress_bar(length):
    return ProgressBar(widgets=[Percentage(),
                                ' (', Counter(), ' of %s) ' % length,
                                Bar()],
                       max_value=length).start()
        
        
def save_classified_reads(name, reads, out_directory):
    with open(os.path.join(out_directory, name), 'w+') as f:
        for read in reads:
            f.write('%s\n' % read)


def load_mapped_reads(name):
    """
    Reads flat text files with unaligned read names and puts them into a dictionary
    organized by (lane, side, tile).

    random_selection is set to a fraction between 0.0 and 1.0 when not all reads should be loaded.
    This is used in cases where phiX is not available for alignment and a subset of all reads should be used.

    """
    read_data = defaultdict(list)
    with open('%s' % os.path.join('mapped_reads', name)) as f:
        for line in f:
            record = FastqName(name=line.strip())
            fastq_read = FastqRead(record)
            read_data[fastq_read.region].append(fastq_read)
    return read_data
