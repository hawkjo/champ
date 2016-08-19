import logging
import os

from champ import align, initialize, error, projectinfo, chip, fastqimagealigner
from champ.config import PathInfo

log = logging.getLogger(__name__)


def main(clargs):
    # TODO: Check if preprocessing is done, if not, run the preprocessing command
    # We know which channel phix is in from the YAML file
    metadata = initialize.load(clargs.image_directory)
    h5_filenames = list(filter(lambda x: x.endswith('.h5'), os.listdir(clargs.image_directory)))
    h5_filenames = [os.path.join(clargs.image_directory, filename) for filename in h5_filenames]
    path_info = PathInfo(clargs.image_directory, metadata['mapped_reads'], clargs.perfect_target_name)

    if len(h5_filenames) == 0:
        error.fail("There were no HDF5 files to process. "
                   "Either they just don't exist, or you didn't provide the correct path.")

    # Ensure we have the directories where output will be written
    align.make_output_directories(h5_filenames, path_info)

    log.debug("Loading tile data.")
    sequencing_chip = chip.load(metadata['chip_type'])(metadata['ports_on_right'])

    alignment_tile_data = align.load_read_names(path_info.aligning_read_names_filepath)
    unclassified_tile_data = align.load_read_names(path_info.all_read_names_filepath)
    perfect_tile_data = align.load_read_names(path_info.perfect_read_names)
    all_tile_data = {key: list(set(alignment_tile_data.get(key, []) + unclassified_tile_data.get(key, [])))
                     for key in list(unclassified_tile_data.keys()) + list(alignment_tile_data.keys())}
    log.debug("Tile data loaded.")

    # We use one process per concentration. We could theoretically speed this up since our machine
    # has significantly more cores than the typical number of concentration points, but since it
    # usually finds a result in the first image or two, it's not going to deliver any practical benefits
    log.debug("Loading FastQImageAligner")
    fia = fastqimagealigner.FastqImageAligner()
    fia.load_reads(alignment_tile_data)
    log.debug("Loaded %s points" % sum([len(v) for v in alignment_tile_data.values()]))
    log.debug("FastQImageAligner loaded.")

    if 'end_tiles' not in metadata:
        end_tiles = align.get_end_tiles(h5_filenames, metadata['alignment_channel'], clargs.snr, metadata, sequencing_chip, fia)
        metadata['end_tiles'] = end_tiles
        initialize.update(clargs.image_directory, metadata)
    else:
        log.debug("End tiles already calculated.")
        end_tiles = metadata['end_tiles']

    if not metadata['phix_aligned']:
        align.run(h5_filenames, path_info, clargs.snr, clargs.min_hits, fia, end_tiles, metadata['alignment_channel'],
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

        # Align all protein reads to the protein image
        channel_combo = channel_name + "_unclassified"
        combo_align(h5_filenames, channel_combo, channel_name, path_info, unclassified_tile_data, all_tile_data, metadata, clargs)

        # Align just perfect protein reads to the protein image (less likely, but might be higher quality alignment!)
        # channel_combo = channel_name + "_perfect"
        # combo_align(h5_filenames, channel_combo, channel_name, path_info, perfect_tile_data, all_tile_data, metadata, clargs)


def combo_align(h5_filenames, channel_combo, channel_name, path_info, alignment_tile_data, all_tile_data, metadata, clargs):
    log.info("Aligning %s" % channel_combo)
    if channel_combo not in metadata['protein_channels_aligned']:
        align.run_data_channel(h5_filenames, channel_name, path_info, alignment_tile_data, all_tile_data, metadata, clargs)
        metadata['protein_channels_aligned'].append(channel_combo)
        initialize.update(clargs.image_directory, metadata)
