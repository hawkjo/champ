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
    nd2_filenames = list(filter(lambda x: x.endswith('.nd2'), os.listdir(os.getcwd())))
    experiment = Experiment(clargs.project_name)
    objective = 60
    alignment_parameters = AlignmentParameters(clargs)
    all_tile_data = reads.get_read_names(alignment_parameters.all_read_names_filepath)
    phix_tile_data = reads.get_read_names(os.path.join(experiment.project_name,
                                                       alignment_parameters.aligning_read_names_filepath))
    processes = min(len(nd2_filenames), multiprocessing.cpu_count())
    log.debug("Using %d processes for alignment" % processes)
    pool = multiprocessing.Pool(processes=processes)
    alignment_function = functools.partial(align.run,
                                           alignment_parameters,
                                           phix_tile_data,
                                           all_tile_data,
                                           experiment,
                                           objective)
    pool.map_async(alignment_function, nd2_filenames).get(sys.maxint)
