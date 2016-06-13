from chimp.config import AlignmentParameters, Experiment
import os
import logging
from chimp import align

log = logging.getLogger(__name__)


def main(clargs):
    h5_filenames = list(filter(lambda x: x.endswith('.h5'), os.listdir(clargs.image_directory)))
    h5_filenames = [os.path.join(clargs.image_directory, filename) for filename in h5_filenames]
    experiment = Experiment(clargs.project_name)
    alignment_parameters = AlignmentParameters(clargs)
    log.debug("Loading tile data.")
    alignment_tile_data = align.load_read_names(alignment_parameters.aligning_read_names_filepath)
    unclassified_tile_data = align.load_read_names(alignment_parameters.all_read_names_filepath)
    all_tile_data = {key: list(set(alignment_tile_data.get(key, []) + unclassified_tile_data.get(key, [])))
                     for key in list(unclassified_tile_data.keys()) + list(alignment_tile_data.keys())}
    log.debug("Tile data loaded.")

    if not clargs.second_channel:
        align.run(h5_filenames, alignment_parameters, alignment_tile_data, all_tile_data, experiment, clargs)
    else:
        align.run_second_channel(h5_filenames, alignment_parameters, all_tile_data, experiment, clargs)
