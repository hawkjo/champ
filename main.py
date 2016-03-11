"""
Chip-Hybridized Interaction Mapping Platform

Usage:
  chimp preprocess [-v | -vv | -vvv]
  chimp align [--chip_id] [--objective] [--min_hits] [--min_tile] [--max_tile] [--fq_w_estimate] [--rotation_estimate] [--snr_threshold] [--index_offset] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

"""
import matplotlib
matplotlib.use('agg')
from docopt import docopt
from controller.align import align
from controller.preprocess import fitsify
import logging


def main(args=None):
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

    if arguments['preprocess']:
        fitsify(arguments)

    if arguments['align']:
        align(arguments)


if __name__ == '__main__':
    main()
