from chimp.config import AlignmentParameters, Experiment
from error import quit
from chimp.process_nd2_im import process_fig


def align(clargs):
    alignment_parameters = AlignmentParameters(clargs)
    file_structure = Experiment()
    process_fig(alignment_parameters, 'fast', clargs.)