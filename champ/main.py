"""
Chip-Hybridized Affinity Mapping Platform

Usage:
  champ map FASTQ_DIRECTORY OUTPUT_DIRECTORY PATHS_TO_BAMFILES ... [--force] [-v | -vv | -vvv]
  champ init IMAGE_DIRECTORY CHIP_NAME READ_NAMES_DIRECTORY ALIGNMENT_CHANNEL LDA_WEIGHTS [--perfect-target-name=PERFECT_TARGET_NAME] [--alternate-perfect-reads=ALTERNATE_PERFECT_READS] [--alternate-good-reads=ALTERNATE_GOOD_READS] [--alternate-fiducial-reads=ALTERNATE_FIDUCIAL_READS] [--microns-per-pixel=0.266666666] [--chip=miseq] [--ports-on-right] [--flipud] [--fliplr] [-v | -vv | -vvv ]
  champ align IMAGE_DIRECTORY [--rotation-adjustment] [--min-hits=MIN_HITS] [--snr=SNR] [--make-pdfs] [--fiducial-only] [-v | -vv | -vvv]
  champ info IMAGE_DIRECTORY

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  map           Maps all the reads in the fastq files. This needs to be done before any other processing
  init          Stores some metadata about a particular experiment
  align         Determines the sequence of fluorescent points in the microscope data. Preprocesses images if not already done.
  info          View the metadata associated with an experiment

"""
import logging
import os
from champ.config import CommandLineArguments
from champ.constants import VERSION
from champ.controller import align, initialize, mapreads, info
from docopt import docopt


def main(**kwargs):
    docopt_args = docopt(__doc__, version=VERSION)
    arguments = CommandLineArguments(docopt_args, os.getcwd())

    log = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(arguments.log_level)
    log.debug(docopt_args)

    commands = {'align': align,
                'init': initialize,
                'map': mapreads,
                'info': info}

    commands[arguments.command].main(arguments)


if __name__ == '__main__':
    main()
