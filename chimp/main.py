"""
Chip-Hybridized Interaction Mapping Platform

Usage:
  chimp bowtie PATH_TO_FASTA [-v | -vv | -vvv]
  chimp readsort FASTQ_DIRECTORY PATHS_TO_BAMFILES ... [-o sorted_files_directory] [-v | -vv | -vvv]
  chimp align [--min_hits] [--tile_width_estimate] [--rotation_estimate] [--snr_threshold] [--index_offset] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  test      ensures that all dependencies are installed
  bowtie    convenience wrapper for a bowtie2 command to create files necessary to classify reads in fastq files
  readsort  classifies all the reads in the fastq files, typically for separating phiX from your reads
  align     maps reads from the high-throughput sequencer to fluorescent points in microscope image data

"""
from controller import align, bowtie, readsort
from docopt import docopt
import error
import logging


def main(**kwargs):
    arguments = docopt(__doc__, version='0.0.1')

    # configure the logger
    log = logging.getLogger()
    log.addHandler(logging.StreamHandler())
    log_level = {0: logging.FATAL,
                 1: logging.WARN,
                 2: logging.INFO,
                 3: logging.DEBUG}
    # default to silent if the user supplies no verbosity setting
    # if they give us an invalid setting (four or more v's) set to debug mode
    log.setLevel(log_level.get(arguments.get('-v', 0), 3))

    # make some space to distinguish log messages from command prompt
    for _ in range(3):
        log.info('')

    try:

        if arguments['readsort']:
            readsort.main(arguments)
        if arguments['bowtie']:
            raise NotImplementedError("We haven't set up the bowtie thing yet. Sorry.")
        if arguments['align']:
            align.main(arguments)

    except Exception as e:
        error.fail(str(e))


if __name__ == '__main__':
    main()
