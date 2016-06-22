import logging
from chimp.fastq import FastqFiles
from chimp import fastq
from chimp import error
import os

log = logging.getLogger(__name__)


def main(clargs):
    """
    Parses fastq files and creates text files containing read names that belong to each source
    Typically this is just used to separate phiX reads from everything else.

    """
    # hardcode the output directory for mapped reads
    out_directory = os.path.join(clargs.fastq_directory, 'mapped_reads')
    # validate and/or create directories
    if not os.path.isdir(clargs.fastq_directory):
        error.fail("The given fastq directory does not exist.")
    if not os.path.isdir(out_directory):
        os.makedirs(out_directory)
    filenames = [os.path.join(clargs.fastq_directory, filename)
                 for filename in os.listdir(clargs.fastq_directory)]

    fastq_files = FastqFiles(filenames)
    all_classified_reads = set()

    if not clargs.force:
        if not fastq.safe_to_classify(clargs.bamfiles, out_directory):
            error.fail("Some or all reads have already been mapped. If you want to force the read mapping "
                       "to be redone, rerun the last command with --force")

    classified_reads = fastq.classify_all_reads(clargs.bamfiles, fastq_files)
    for name, reads in classified_reads.items():
        # write the reads to disk for later use
        fastq.save_classified_reads(name, reads, out_directory)
        # keep a single set with all the reads so that we can figure out which reads aren't classified
        all_classified_reads.update(reads)

    # now figure out which reads remain unclassified and save them to a catch-all bucket
    log.info('Finding reads that were not classified.')
    all_unclassified_reads = fastq.load_unclassified_reads(fastq_files, all_classified_reads)
    fastq.save_classified_reads('unclassified', all_unclassified_reads, out_directory)
