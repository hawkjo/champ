from champ.config import AlignmentParameters, Experiment
import os
import logging
from champ import align, initialize

log = logging.getLogger(__name__)


def main(clargs):
    # TODO: KeyError: 'SECOND_CHANNEL_NAME'
    # TODO: Check if preprocessing is done, if not, run the preprocessing command
    # TODO: for each channel, determine if alignment is complete, and if not, align that channel, starting with phix first
    # We know which channel phix is in from the YAML file
    # TODO: if --phix-only, don't proceed with 2nd channels
    # TODO: add auto-elbow-grease, a technique to align images with an abnormally low number of clusters
    metadata = initialize.load(clargs.image_directory)
    h5_filenames = list(filter(lambda x: x.endswith('.h5'), os.listdir(clargs.image_directory)))
    h5_filenames = [os.path.join(clargs.image_directory, filename) for filename in h5_filenames]
    experiment = Experiment(clargs.image_directory)
    alignment_parameters = AlignmentParameters(clargs, metadata['mapped_reads'])
    log.debug("Loading tile data.")
    alignment_tile_data = align.load_read_names(alignment_parameters.aligning_read_names_filepath)
    unclassified_tile_data = align.load_read_names(alignment_parameters.all_read_names_filepath)
    all_tile_data = {key: list(set(alignment_tile_data.get(key, []) + unclassified_tile_data.get(key, [])))
                     for key in list(unclassified_tile_data.keys()) + list(alignment_tile_data.keys())}
    log.debug("Tile data loaded.")

    if not clargs.second_channel:
        align.run(h5_filenames, alignment_parameters, alignment_tile_data, all_tile_data, experiment, metadata)
    else:
        align.run_second_channel(h5_filenames, alignment_parameters, all_tile_data, experiment, clargs)
