import logging

import os
from champ import align, initialize
from champ import projectinfo
from champ.config import AlignmentParameters, Experiment

log = logging.getLogger(__name__)


def main(clargs):
    # TODO: Check if preprocessing is done, if not, run the preprocessing command
    # TODO: for each channel, determine if alignment is complete, and if not, align that channel, starting with phix
    # We know which channel phix is in from the YAML file
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

    align.run(h5_filenames, alignment_parameters, alignment_tile_data, all_tile_data, experiment, metadata, clargs.make_pdfs)
    if not clargs.phix_only:
        protein_channels = [channel for channel in projectinfo.load_channels(clargs.image_directory)
                            if channel != metadata['alignment_channel']]
        log.debug("Protein channels found: %s" % ", ".join(protein_channels))
        for channel_name in protein_channels:
            log.debug("Aligning protein channel: %s" % channel_name)
            align.run_data_channel(h5_filenames, channel_name, alignment_parameters, alignment_tile_data,
                                   all_tile_data, experiment, metadata, clargs)
