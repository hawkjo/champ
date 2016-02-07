"""Chip-Hybridized Interaction Mapping Platform

Usage:
  chimp <command> [-v | -vv | -vvv]

Commands:
  align         Creates alignments from raw images and NGS sequence data.

Options:
  -h --help     Show this screen.
  --version     Show version.

"""
from docopt import docopt
from error import quit
from preprocess.main import run as align
import logging
import matplotlib
matplotlib.use('agg')


if __name__ == '__main__':
    arguments = docopt(__doc__, version='0.0.1')

    # configure the logger
    log = logging.getLogger()
    log.addHandler(logging.StreamHandler())
    log_level = {0: logging.FATAL,
                 1: logging.WARN,
                 2: logging.INFO,
                 3: logging.DEBUG}
    log.setLevel(log_level.get(arguments.get('-v', 0), 3))

    # parse the command
    commands = {'align': align}
    command_name = arguments['<command>']
    if command_name not in commands:
        quit("Invalid command.")

    # run the command
    commands[command_name]()
