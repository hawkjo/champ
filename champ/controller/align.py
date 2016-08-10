import logging
import os

from champ import align, initialize, error, projectinfo, chip, fastqimagealigner
from champ.config import OutputParameters

log = logging.getLogger(__name__)


def main(clargs):
    # TODO: Check if preprocessing is done, if not, run the preprocessing command
    # We know which channel phix is in from the YAML file
    metadata = initialize.load(clargs.image_directory)
    h5_filenames = list(filter(lambda x: x.endswith('.h5'), os.listdir(clargs.image_directory)))
    h5_filenames = [os.path.join(clargs.image_directory, filename) for filename in h5_filenames]
    output_parameters = OutputParameters(clargs.image_directory, metadata['mapped_reads'])

    if len(h5_filenames) == 0:
        error.fail("There were no HDF5 files to process. "
                   "Either they just don't exist, or you didn't provide the correct path.")

    # Ensure we have the directories where output will be written
    align.make_output_directories(h5_filenames, output_parameters)

    log.debug("Loading tile data.")
    sequencing_chip = chip.load(metadata['chip_type'])(metadata['ports_on_right'])
    alignment_tile_data = align.load_read_names(output_parameters.aligning_read_names_filepath)
    unclassified_tile_data = align.load_read_names(output_parameters.all_read_names_filepath)
    all_tile_data = {key: list(set(alignment_tile_data.get(key, []) + unclassified_tile_data.get(key, [])))
                     for key in list(unclassified_tile_data.keys()) + list(alignment_tile_data.keys())}
    log.debug("Tile data loaded.")

    # We use one process per concentration. We could theoretically speed this up since our machine
    # has significantly more cores than the typical number of concentration points, but since it
    # usually finds a result in the first image or two, it's not going to deliver any practical benefits
    fia = fastqimagealigner.FastqImageAligner()
    fia.load_reads(alignment_tile_data)

    if 'end_tiles' not in metadata:
        end_tiles = align.get_end_tiles(h5_filenames, metadata['alignment_channel'], clargs.snr, metadata, sequencing_chip, fia)
        metadata['end_tiles'] = end_tiles
        initialize.update(clargs.image_directory, metadata)
    else:
        log.debug("End tiles already calculated.")
        end_tiles = metadata['end_tiles']

    if not metadata['phix_aligned']:
        align.run(h5_filenames, output_parameters, clargs.snr, clargs.min_hits, fia, end_tiles, metadata['alignment_channel'],
                  all_tile_data, metadata, clargs.make_pdfs, sequencing_chip)
        metadata['phix_aligned'] = True
        initialize.update(clargs.image_directory, metadata)
    else:
        log.debug("Phix already aligned.")

    protein_channels = [channel for channel in projectinfo.load_channels(clargs.image_directory) if channel != metadata['alignment_channel']]
    if protein_channels:
        log.debug("Protein channels found: %s" % ", ".join(protein_channels))
    else:
        # protein is in phix channel, hopefully?
        log.warn("No protein channels detected. Assuming protein is in phiX channel: %s" % [metadata['alignment_channel']])
        protein_channels = [metadata['alignment_channel']]
    for channel_name in protein_channels:
        if channel_name not in metadata['protein_channels_aligned']:
            log.debug("Aligning protein channel: %s" % channel_name)
            align.run_data_channel(h5_filenames, channel_name, output_parameters, alignment_tile_data,
                                   all_tile_data, metadata, clargs)
            metadata['protein_channels_aligned'].append(channel_name)
            initialize.update(clargs.image_directory, metadata)
