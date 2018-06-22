import logging
from champ import genome, error
import os

log = logging.getLogger(__name__)


def main(clargs):
    if not os.path.isdir(clargs.fastq_directory):
        error.fail("The given fastq directory does not exist.")
    genome.build_genomic_bamfile(clargs.fastq_directory, clargs.output_filename)
