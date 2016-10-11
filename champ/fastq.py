import gzip
import logging
import os


log = logging.getLogger(__name__)



def classify_fastq_reads(classifier, fastq_files):
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
        classifier = FastqReadClassifier(path)
        name, reads = classify_fastq_reads(classifier, fastq_files)
        assert name != 'unclassified', '"unclassified" cannot be used as a fastq read classifier name'
        classified_reads[name] = reads
    return classified_reads


def safe_to_classify(classifier_paths, out_directory):
    # Will we overwrite existing classifications? If so, we're probably doing redundant work unless
    # the classification process was halted before it finished
    for path in classifier_paths:
        classifier = FastqReadClassifier(path)
        filename = os.path.join(out_directory, classifier.name)
        if os.path.isfile(filename):
            return False
    return True


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
