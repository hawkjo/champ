from config import AlignmentParameters
from error import quit
import os
from process_nd2_im import process_fig


def guess_chip_id(directory):
    """ Tries to figure out which directory in the project directory contains the sequencing data. """
    for name in os.listdir(directory):
        possible_chip_dir = os.path.sep.join((directory, name))
        if os.path.isdir(possible_chip_dir):
            files = [f for f in os.listdir(possible_chip_dir)]
            if 'all_fastqs' in files:
                return name.rstrip(os.sep)
    return None
