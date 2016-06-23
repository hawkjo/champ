"""
Chip-Hybridized Interaction Mapping Platform

Usage:
  chimp convert TIF_DIRECTORIES ... [--flipud] [--fliplr] [-v | -vv | -vvv]
  chimp preprocess IMAGE_DIRECTORY [-v | -vv | -vvv ]
  chimp map FASTQ_DIRECTORY OUTPUT_DIRECTORY PATHS_TO_BAMFILES ... [--force] [-v | -vv | -vvv]
  chimp align ALIGNMENT_CHANNEL IMAGE_DIRECTORY PROJECT_NAME MICRONS_PER_PIXEL [--chip=miseq] [--second-channel SECOND_CHANNEL_NAME] [--ports-on-right] [--min-hits MIN_HITS] [--snr-threshold SNR] [--make-pdfs] [-v | -vv | -vvv]
  chimp intensity [-v | -vv | -vvv]
  chimp info IMAGE_DIRECTORY [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  convert       creates an HDF5-formatted file from OME-TIFF files
  map           maps all the reads in the fastq files, typically for separating phiX
  preprocess    defines where points are in the microscope image data
  align         maps reads from the high-throughput sequencer to fluorescent
                points in microscope image data
  intensity     determines boundaries of clusters and assigns intensities to sequences
  info          brief summary of the data

"""
from chimp.controller import align, preprocess, mapreads, convert, intensity, info
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
    log.debug(str(arguments.__dict__))

    # make some space to distinguish log messages from command prompt
    for _ in range(2):
        log.info('')

    commands = {'align': align,
                'preprocess': preprocess,
                'map': mapreads,
                'convert': convert,
                'intensity': intensity,
                'info': info}

    commands[arguments.command].main(arguments)


if __name__ == '__main__':
    main()
