"""
Chip-Hybridized Affinity Mapping Platform

Usage:
  champ map FASTQ_DIRECTORY OUTPUT_DIRECTORY [--log-p-file=LOG_P_FILE] [--target-sequence-file=TARGET_SEQUENCE_FILE] [--phix-bowtie=PHIX_BOWTIE] [--min-len=MIN_LEN] [--max-len=MAX_LEN] [--include-side-1] [-v | -vv | -vvv]
  champ init IMAGE_DIRECTORY READ_NAMES_DIRECTORY [ALIGNMENT_CHANNEL] [--perfect-target-name=PERFECT_TARGET_NAME] [--alternate-perfect-reads=ALTERNATE_PERFECT_READS] [--alternate-good-reads=ALTERNATE_GOOD_READS] [--alternate-fiducial-reads=ALTERNATE_FIDUCIAL_READS] [--microns-per-pixel=0.266666666] [--chip=miseq] [--ports-on-right] [--flipud] [--fliplr] [-v | -vv | -vvv ]
  champ h5 IMAGE_DIRECTORY [-v | -vv | -vvv]
  champ align IMAGE_DIRECTORY [--rotation-adjustment=ROTATION_ADJUSTMENT] [--min-hits=MIN_HITS] [--snr=SNR] [--process-limit=PROCESS_LIMIT] [--make-pdfs] [--fiducial-only] [-v | -vv | -vvv]
  champ info IMAGE_DIRECTORY
  champ notebooks

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  map           Maps all the reads in the fastq files. This needs to be done before any other processing
  init          Stores some metadata about a particular experiment
  h5            Looks into all directories in IMAGE_DIRECTORY and converts any valid TIFFs it finds into HDF5 files
  align         Determines the sequence of fluorescent points in the microscope data. Preprocesses images if not already done
  info          View the metadata associated with an experiment
  notebooks     Creates copies of the standard Jupyter analysis notebooks in the current directory

"""
import logging
import os
from champ.config import CommandLineArguments
from champ.constants import VERSION
from champ.controller import align, initialize, h5, mapreads, info, notebooks
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
                'h5': h5,
                'map': mapreads,
                'info': info,
                'notebooks': notebooks}

    commands[arguments.command].main(arguments)


if __name__ == '__main__':
    main()
