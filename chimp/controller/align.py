from chimp.config import AlignmentParameters, Experiment
import os
import logging
import sys
from chimp import align
import multiprocessing
import functools

log = logging.getLogger(__name__)


def main(clargs):
    h5_filenames = list(filter(lambda x: x.endswith('.h5'), os.listdir(os.getcwd())))
    experiment = Experiment(clargs.project_name)
    um_per_pixel = 0.26666666
    alignment_parameters = AlignmentParameters(clargs)
    all_tile_data = align.load_read_names(alignment_parameters.all_read_names_filepath)
    phix_tile_data = align.load_read_names(alignment_parameters.aligning_read_names_filepath)

    processes = min(len(h5_filenames), multiprocessing.cpu_count())
    log.debug("Using %d processes for alignment" % processes)
    pool = multiprocessing.Pool(processes=processes)
    alignment_function = functools.partial(align.run,
                                           alignment_parameters,
                                           phix_tile_data,
                                           all_tile_data,
                                           experiment,
                                           um_per_pixel,
                                           clargs.alignment_channel)

    pool.map_async(alignment_function, h5_filenames).get(sys.maxint)
