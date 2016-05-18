"""
Chip-Hybridized Interaction Mapping Platform

Usage:
  chimp convert TIF_DIRECTORIES ... [--flipud] [--fliplr] [-v | -vv | -vvv]
  chimp preprocess IMAGE_DIRECTORY [-v | -vv | -vvv ]
  chimp map FASTQ_DIRECTORY PATHS_TO_BAMFILES ... [-v | -vv | -vvv]
  chimp align ALIGNMENT_CHANNEL IMAGE_DIRECTORY PROJECT_NAME [--min-hits] [--snr-threshold] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  convert       creates an HDF5-formatted file from OME-TIFF files
  map           maps all the reads in the fastq files, typically for separating phiX
  preprocess    defines where points are in the microscope image data
  align         maps reads from the high-throughput sequencer to fluorescent
                points in microscope image data

"""
from chimp.controller import align, preprocess, mapreads, convert
from docopt import docopt
import logging
from chimp.config import CommandLineArguments
from chimp.constants import VERSION
import os


def main(**kwargs):
    arguments = CommandLineArguments(docopt(__doc__, version=VERSION), os.getcwd())

    log = logging.getLogger()
    log.addHandler(logging.StreamHandler())
    log.setLevel(arguments.log_level)

    # make some space to distinguish log messages from command prompt
    for _ in range(2):
        log.info('')

    commands = {'align': align.main,
                'second': align.second,
                'preprocess': preprocess.main,
                'map': mapreads.main,
                'convert': convert.main
                }
    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
