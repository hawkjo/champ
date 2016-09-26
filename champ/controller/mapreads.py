import logging
from champ.fastq import FastqFiles
from champ import fastq
from champ import error
import os

log = logging.getLogger(__name__)


def main(clargs):
    """
    Parses fastq files and creates text files containing read names, with one file per source
    For example, all reads that are part of the phiX genome will go into a file called "phix"
    These sources are provided by the user, in the form of a list of BAM files.

    """
    # validate and/or create directories
    if not os.path.isdir(clargs.fastq_directory):
        error.fail("The given fastq directory does not exist.")
    if not os.path.isdir(clargs.output_directory):
        os.makedirs(clargs.output_directory)
    filenames = [os.path.join(clargs.fastq_directory, filename)
                 for filename in os.listdir(clargs.fastq_directory)]

    fastq_files = FastqFiles(filenames)
    all_classified_reads = set()

    if not clargs.force and not fastq.safe_to_classify(clargs.bamfiles, clargs.output_directory):
            error.fail("Some or all reads have already been mapped. If you want to force the read mapping "
                       "to be redone, rerun the last command with --force")

    # create all_read_names.txt

    classified_reads = fastq.classify_all_reads(clargs.bamfiles, fastq_files)
    for name, reads in classified_reads.items():
        # write the reads to disk for later use
        fastq.save_classified_reads(name, reads, clargs.output_directory)
        # keep a single set with all the reads so that we can figure out which reads aren't classified
        all_classified_reads.update(reads)

    # do that ML thing to make read_names_by_sequence.txt
    # create perfect_target_read_names
    # create target_read_names
