"""
Chip-Hybridized Interaction Mapping Platform

Usage:
  chimp bowtie PATH_TO_FASTA [-v | -vv | -vvv]
  chimp sort FASTQ_DIRECTORY PATHS_TO_BAMFILES ... [-o sorted_files_directory] [-v | -vv | -vvv]
  chimp align [--min_hits] [--tile_width_estimate] [--rotation_estimate] [--snr_threshold] [--index_offset] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  bowtie    convenience wrapper for a bowtie2 command to create files necessary to classify reads in fastq files
  sort      classifies all the reads in the fastq files, typically for separating phiX from your reads
  align     maps reads from the high-throughput sequencer to fluorescent points in microscope image data

"""
from chimp.controller import align, bowtie, readsort
from docopt import docopt
from chimp import error
import logging
from chimp.model.clargs import CommandLineArguments
from chimp.model.constants import VERSION
import os


def main(**kwargs):
    arguments = CommandLineArguments(docopt(__doc__, version=VERSION), os.getcwd())

        # configure the logger
    log = logging.getLogger()
    log.addHandler(logging.StreamHandler())
    log.setLevel(arguments.log_level)

    # make some space to distinguish log messages from command prompt
    for _ in range(3):
        log.info('')

    commands = {'align': align.main,
                'bowtie': lambda *x: None,
                'sort': readsort.main,
                }
    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
