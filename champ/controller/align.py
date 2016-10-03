import logging
import os
from champ import align, initialize, error, projectinfo, chip, fastqimagealigner, convert, fits
from champ.config import PathInfo
from champ.fastqtilercs import FastqTileRCs

log = logging.getLogger(__name__)


def preprocess(clargs, metadata):
    log.debug("Preprocessing images.")
    paths = convert.get_all_tif_paths(clargs.image_directory)
    # directories will have ".h5" appended to them to come up with the HDF5 names
    # tifs are relative paths to each tif file
    log.debug("About to convert TIFs to HDF5.")
    convert.main(paths, metadata['flipud'], metadata['fliplr'])
    log.debug("Done converting TIFs to HDF5.")
    log.debug("Fitsifying images from HDF5 files.")
    fits.main(clargs.image_directory)
    metadata['preprocessed'] = True
    initialize.update(clargs.image_directory, metadata)


def load_filenames(image_directory):
    h5_filenames = list(filter(lambda x: x.endswith('.h5'), os.listdir(image_directory)))
    return [os.path.join(image_directory, filename) for filename in h5_filenames]


def load_fastq_tiles(tile_data, microns_per_pixel):
    fastq_tiles = {}
    for tile_key, read_names in tile_data.items():
        fastq_tiles[tile_key] = FastqTileRCs(tile_key, read_names, microns_per_pixel)
    return fastq_tiles


def main(clargs):
    metadata = initialize.load(clargs.image_directory)
    if 'preprocessed' not in metadata or not metadata['preprocessed']:
        for filename in load_filenames(clargs.image_directory):
            log.warn("Deleting (probably invalid) existing HDF5 file and recreating it: %s" % filename)
            os.unlink(filename)
        preprocess(clargs, metadata)

    h5_filenames = load_filenames(clargs.image_directory)
    if len(h5_filenames) == 0:
        error.fail("There were no HDF5 files to process. You must have deleted or moved them after preprocessing them.")

    path_info = PathInfo(clargs.image_directory,
                         metadata['mapped_reads'],
                         metadata['perfect_target_name'],
                         metadata['alternate_fiducial_reads'],
                         metadata['alternate_perfect_target_reads_filename'],
                         metadata['alternate_good_target_reads_filename'])
    # Ensure we have the directories where output will be written
    align.make_output_directories(h5_filenames, path_info)

    log.debug("Loading tile data.")
    sequencing_chip = chip.load(metadata['chip_type'])(metadata['ports_on_right'])

    alignment_tile_data = align.load_read_names(path_info.aligning_read_names_filepath)
    all_tile_data = align.load_read_names(path_info.all_read_names_filepath)
    log.debug("Tile data loaded.")

    # We use one process per concentration. We could theoretically speed this up since our machine
    # has significantly more cores than the typical number of concentration points, but since it
    # usually finds a result in the first image or two, it's not going to deliver any practical benefits
    log.debug("Loading FastQImageAligner")
    fastq_tiles = load_fastq_tiles(alignment_tile_data, metadata['microns_per_pixel'])
    log.debug("Loaded %s points" % sum([len(v) for v in alignment_tile_data.values()]))
    log.debug("FastQImageAligner loaded.")

    if 'end_tiles' not in metadata:
        # WARNING: Last parameter fastq_tiles is not correct or expected!!!
        end_tiles = align.get_end_tiles(h5_filenames, metadata['alignment_channel'], clargs.snr, metadata, sequencing_chip, fastq_tiles)
        metadata['end_tiles'] = end_tiles
        initialize.update(clargs.image_directory, metadata)
    else:
        log.debug("End tiles already calculated.")
        end_tiles = metadata['end_tiles']

    if not metadata['phix_aligned']:
        align.align_fiducial(h5_filenames, path_info, clargs.snr, clargs.min_hits, fastq_tiles, end_tiles,
                             metadata['alignment_channel'], all_tile_data, metadata, clargs.make_pdfs, sequencing_chip)
        metadata['phix_aligned'] = True
        initialize.update(clargs.image_directory, metadata)
    else:
        log.debug("Phix already aligned.")

    if clargs.fiducial_only:
        # the user doesn't want us to align the protein channels
        exit(0)

    protein_channels = [channel for channel in projectinfo.load_channels(clargs.image_directory) if channel != metadata['alignment_channel']]
    if protein_channels:
        log.debug("Protein channels found: %s" % ", ".join(protein_channels))
    else:
        # protein is in phix channel, hopefully?
        log.warn("No protein channels detected. Assuming protein is in phiX channel: %s" % [metadata['alignment_channel']])
        protein_channels = [metadata['alignment_channel']]

    log.debug("Loading target read data")
    perfect_tile_data = align.load_read_names(path_info.perfect_read_names)
    on_target_tile_data = align.load_read_names(path_info.on_target_read_names)
    log.debug("Target read data loaded.")

    for channel_name in protein_channels:
        # Attempt to precision align protein channels using the phix channel alignment as a starting point.
        # Not all experiments have "on target" or "perfect target" reads - that only applies to CRISPR systems (at the time of this writing anyway)

        if on_target_tile_data:
            channel_combo = channel_name + "_on_target"
            combo_align(h5_filenames, channel_combo, channel_name, path_info, on_target_tile_data, all_tile_data, metadata, clargs)

        if perfect_tile_data:
            channel_combo = channel_name + "_perfect_target"
            combo_align(h5_filenames, channel_combo, channel_name, path_info, perfect_tile_data, all_tile_data, metadata, clargs)


def combo_align(h5_filenames, channel_combo, channel_name, path_info, alignment_tile_data, all_tile_data, metadata, clargs):
    log.info("Aligning %s" % channel_combo)
    if channel_combo not in metadata['protein_channels_aligned']:
        align.run_data_channel(h5_filenames, channel_name, path_info, alignment_tile_data, all_tile_data, metadata, clargs)
        metadata['protein_channels_aligned'].append(channel_combo)
        initialize.update(clargs.image_directory, metadata)
