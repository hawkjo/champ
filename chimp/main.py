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
from chimp.controller import align, bowtie, readsort
from docopt import docopt
import error
import logging
from chimp.model.clargs import CommandLineArguments
import os


def main(**kwargs):
    arguments = CommandLineArguments(docopt(__doc__, version='0.0.1'), os.getcwd())

        # configure the logger
    log = logging.getLogger()
    log.addHandler(logging.StreamHandler())
    log.setLevel(arguments.log_level)

    # make some space to distinguish log messages from command prompt
    for _ in range(3):
        log.info('')

    commands = {'align': align.main,
                'bowtie': lambda *x: None,
                'readsort': readsort.main,
                }
    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
