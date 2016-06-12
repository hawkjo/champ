from Bio import SeqIO
from collections import namedtuple
import gzip
import logging
import os
import pysam
import subprocess


log = logging.getLogger(__name__)


FastqName = namedtuple('FastqRead', ['name'])


class FastqAlignmentRead(object):
    """ Wraps the raw data about a single DNA read that we receive from Illumina.
        Discards the sequence and quality data to conserve memory. """
    __slots__ = ('name', '_lane', '_side', '_tile')

    def __init__(self, record):
        self.name = record.name
        self._lane = None
        self._side = None
        self._tile = None

    @property
    def region(self):
        return self.lane, self.side, self.tile

    @property
    def lane(self):
        return int(self._lookup_name_data(4))

    @property
    def tile(self):
        return int(self._lookup_name_data(3)[-2:])

    @property
    def side(self):
        return int(self._lookup_name_data(3)[0])

    @property
    def row(self):
        return int(self._lookup_name_data(1))

    @property
    def column(self):
        return int(self._lookup_name_data(2))

    def _lookup_name_data(self, index):
        return self.name.rsplit(':')[-index]


class FastqFiles(object):
    """ Sorts compressed FastQ files provided to us from the Illumina sequencer. """
    def __init__(self, filenames):
        self._filenames = list(self._filter_names(filenames))

    def __iter__(self):
        for f in self._filenames:
            yield f

    def __len__(self):
        return len(self._filenames)

    @property
    def alignment_length(self):
        paired_length = len([(f1, f2) for f1, f2 in self.paired])
        single_length = len([f for f in self.single])
        return paired_length + single_length

    @property
    def paired(self):
        for f1, f2 in self._sort_filenames(paired=True):
            yield f1, f2

    @property
    def single(self):
        for f in self._sort_filenames(paired=False):
            yield f

    def _filter_names(self, data):
        # eliminate filenames that can't possibly be fastq files of interest
        for filename in reversed(data):
            if not filename.endswith('fastq.gz'):
                continue
            if '_I1_' in filename or '_I2_' in filename or '_I1.' in filename or '_I2.' in filename:
                continue
            yield filename

    def _sort_filenames(self, paired=True):
        # yield filenames that are the given type (single or paired)
        for filename in self._filenames:
            if '_R1_' in filename or '_R1.' in filename:
                pair = filename.replace('_R1_', '_R2_').replace('_R1.', '_R2.')
                if paired and pair in self._filenames:
                    yield filename, pair
                elif not paired and pair not in self._filenames:
                    yield filename


def parse_fastq_lines(fh):
    for dna in SeqIO.parse(fh, 'fastq'):
        yield FastqAlignmentRead(dna)


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
        # TODO: ELIMINATE THIS HARD CODED PATH
        with open('/mnt/marble/hdd/home/shared/chimp.log', 'w+') as devnull:
            kwargs = dict(shell=True, stderr=devnull, stdout=devnull)
            subprocess.call(' '.join(command), **kwargs)
            sam_command = 'samtools view -bS /tmp/chimp.sam | samtools sort - /tmp/final'
            subprocess.call(sam_command, **kwargs)
            subprocess.call('samtools index /tmp/final.bam', **kwargs)
            for r in pysam.Samfile('/tmp/final.bam'):
                yield r.qname


def classify_fastq_reads(classifier_path, fastq_files):
    classifier = FastqReadClassifier(classifier_path)
    log.info('Searching reads for %s. This will take a while!' % classifier.name)

    current = 0
    reads = set()
    for file1, file2 in fastq_files.paired:
        for read in classifier.paired_call(file1, file2):
            reads.add(read)
        current += 1
    for file1 in fastq_files.single:
        for read in classifier.single_call(file1):
            reads.add(read)
        current += 1
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
        with gzip.open(filename) as fh:
            yield {r.name for r in parse_fastq_lines(fh)}
        del fh


def load_unclassified_reads(fastq_files, all_classified_reads):
    all_unclassified_reads = set()
    for all_read_names in stream_all_read_names(fastq_files):
        all_unclassified_reads.update(all_read_names.difference(all_classified_reads))
    return all_unclassified_reads

        
def save_classified_reads(name, reads, out_directory):
    with open(os.path.join(out_directory, name), 'w+') as f:
        for read in reads:
            f.write('%s\n' % read)
