"""
Chip-Hybridized Affinity Mapping Platform

Usage:
  champ map FASTQ_DIRECTORY OUTPUT_DIRECTORY PATHS_TO_BAMFILES ... [--force] [-v | -vv | -vvv]
  champ init IMAGE_DIRECTORY CHIP_NAME MAPPED_READS ALIGNMENT_CHANNEL [--microns-per-pixel=0.266666666] [--chip=miseq] [--ports-on-right] [--flipud] [--fliplr] [-v | -vv | -vvv ]
  champ info IMAGE_DIRECTORY [-v | -vv | -vvv]
  champ align IMAGE_DIRECTORY [--phix-only] [--min-hits MIN_HITS] [--snr-threshold SNR] [--make-pdfs] [-v | -vv | -vvv]
  champ kd IMAGE_DIRECTORY TARGET_DATA_FILE TARGET_LABEL OFF_TARGET_LABEL [-v | -vv | -vvv]
  champ preprocess IMAGE_DIRECTORY [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  map           Maps all the reads in the fastq files. This needs to be done before any other processing
  init          Stores some metadata about a particular experiment
  preprocess    Convert TIFs to HDF5 and prepare microscope data for alignment. Only needed for development.
  align         Determines the sequence of fluorescent points in the microscope data. Preprocesses images if not already done.
  kd            Determines boundaries of clusters, assigns intensities to sequences and derives the apparent Kd's
  info          brief summary of the data

"""
from chimp.controller import align, init, mapreads, kd, info, preprocess
from docopt import docopt
import logging
from chimp.config import CommandLineArguments
from chimp.constants import VERSION
import os


def main(**kwargs):
    docopt_args = docopt(__doc__, version=VERSION)
    arguments = CommandLineArguments(docopt_args, os.getcwd())

    log = logging.getLogger()
    log.addHandler(logging.StreamHandler())
    log.setLevel(arguments.log_level)

    # make some space to distinguish log messages from command prompt
    for _ in range(2):
        log.info('')

    commands = {'align': align,
                'preprocess': preprocess,
                'init': init,
                'map': mapreads,
                'kd': kd,
                'info': info}

    commands[arguments.command].main(arguments)


if __name__ == '__main__':
    main()
