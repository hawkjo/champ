from chimp.config import AlignmentParameters, Experiment
import os
import logging
import sys
from chimp import align
import multiprocessing
import functools

log = logging.getLogger(__name__)


def main(clargs):
    h5_filenames = list(filter(lambda x: x.endswith('.h5'), os.listdir(clargs.image_directory)))
    h5_filenames = [os.path.join(clargs.image_directory, filename) for filename in h5_filenames]
    experiment = Experiment(clargs.project_name)
    um_per_pixel = 0.2666666666
    alignment_parameters = AlignmentParameters(clargs)
    log.debug("Loading tile data.")
    phix_tile_data = align.load_read_names(alignment_parameters.aligning_read_names_filepath)
    log.debug("Tile data loaded.")
    processes = min(len(h5_filenames), multiprocessing.cpu_count())
    log.debug("Using %d processes for alignment" % processes)
    align.run(h5_filenames, alignment_parameters, phix_tile_data,
              experiment, um_per_pixel, clargs.alignment_channel)


# def second(clargs):
#     h5_filenames = list(filter(lambda x: x.endswith('.h5'), os.listdir(clargs.image_directory)))
#     h5_filenames = [os.path.join(clargs.image_directory, filename) for filename in h5_filenames]
#     experiment = Experiment(clargs.project_name)
#     um_per_pixel = 0.2666666666
#     alignment_parameters = AlignmentParameters(clargs)
#     log.debug("Loading tile data.")
#     tile_data = align.load_read_names(alignment_parameters.all_read_names_filepath)
#     log.debug("Tile data loaded.")
#     processes = min(len(h5_filenames), multiprocessing.cpu_count())
#     log.debug("Using %d processes for alignment" % processes)
#     pool = multiprocessing.Pool(processes=processes)
#
#     align.process_second_fig(alignment_parameters, base_name, tile_data, um_per_pixel,
#                              experiment, image)
#     alignment_function = functools.partial(align.run,
#                                            alignment_parameters,
#                                            tile_data,
#                                            experiment,
#                                            um_per_pixel,
#                                            clargs.alignment_channel)
#     pool.map_async(alignment_function, h5_filenames).get(sys.maxint)
