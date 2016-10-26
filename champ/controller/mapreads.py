import logging
from champ import error, readmap
import os

log = logging.getLogger(__name__)


def main(clargs):
    # validate and/or create directories
    if not os.path.isdir(clargs.fastq_directory):
        error.fail("The given fastq directory does not exist.")
    if not os.path.isdir(clargs.output_directory):
        os.makedirs(clargs.output_directory)

    readmap.main(clargs)
