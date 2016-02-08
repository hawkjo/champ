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
            if 'all_fastqs' in files and 'phiX_mappings' in files:
                return name.rstrip(os.sep)
    return None


def align(command_line_args):
    cwd = os.getcwd()
    presumptive_chip_id = guess_chip_id(cwd)
    alignment_parameters = AlignmentParameters(command_line_args, cwd, presumptive_chip_id)
    if alignment_parameters.chip_id is None:
        quit("No chip ID was given and the chip ID could not be guessed. You must pass it in manually with --chip_id")

    # Now John needs to explain to Jim how end tile selection works
    process_fig(alignment_parameters, 'fast', '15-11-18_SA15243_Cascade-TA_1nM-007.nd2', 367)
