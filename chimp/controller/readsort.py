import error
import fastq
import logging
from model.fastq import FastqFiles
import os

log = logging.getLogger(__name__)


def main(arguments):
    """
    Parses fastq files and creates text files containing read names that belong to each source
    Typically this is just used to separate phiX reads from everything else.

    """
    # validate and/or create directories
    out_directory = arguments['-o'] if arguments['-o'] else os.path.join(os.getcwd(), 'classified_reads')
    log.info("outdir %s" % out_directory)
    if not os.path.isdir(out_directory):
        os.makedirs(out_directory)
    if not os.path.isdir(arguments['FASTQ_DIRECTORY']):
        error.fail("The given fastq directory does not exist.")
    filenames = [os.path.join(arguments['FASTQ_DIRECTORY'], filename)
                 for filename in os.listdir(arguments['FASTQ_DIRECTORY'])]

    fastq_files = FastqFiles(filenames)
    all_classified_reads = set()

    classified_reads = fastq.classify_all_reads(arguments['PATHS_TO_BAMFILES'], fastq_files)
    for name, reads in classified_reads.items():
        # write the reads to disk for later use
        fastq.save_classified_reads(name, reads, out_directory)
        # keep a single set with all the reads so that we can figure out which reads aren't classified
        all_classified_reads.update(reads)

    # now figure out which reads remain unclassified and save them to a catch-all bucket
    log.info('Finding reads that were not classified.')
    all_unclassified_reads = fastq.load_unclassified_reads(fastq_files, all_classified_reads)
    fastq.save_classified_reads('unclassified', all_unclassified_reads, out_directory)
