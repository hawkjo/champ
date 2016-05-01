from chimp.config import AlignmentParameters, Experiment
from chimp import reads
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
    um_per_pixel = 0.27
    alignment_parameters = AlignmentParameters(clargs)
    all_tile_data = reads.get_read_names(alignment_parameters.all_read_names_filepath)
    phix_tile_data = reads.get_read_names(alignment_parameters.aligning_read_names_filepath)

    # Jim's laptop only
    for h5_filename in h5_filenames:
        align.run(alignment_parameters, phix_tile_data,
                  all_tile_data, experiment, um_per_pixel, h5_filename)

    # Actual code we should use
    # processes = min(len(h5_filenames), multiprocessing.cpu_count())
    # log.debug("Using %d processes for alignment" % processes)
    # pool = multiprocessing.Pool(processes=processes)
    # alignment_function = functools.partial(align.run,
    #                                        alignment_parameters,
    #                                        phix_tile_data,
    #                                        all_tile_data,
    #                                        experiment,
    #                                        objective)
    # pool.map_async(alignment_function, h5_filenames).get(sys.maxint)
