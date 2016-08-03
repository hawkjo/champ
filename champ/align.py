import matplotlib
matplotlib.use('Agg')
from champ.grid import GridImages
from champ import plotting, error, chip
from collections import Counter, defaultdict
import fastqimagealigner
import functools
import h5py
import logging
import multiprocessing
from multiprocessing import Manager
import os
import sys
import re
from copy import deepcopy

log = logging.getLogger(__name__)
stats_regex = re.compile(r'''^(\w+)_(?P<row>\d+)_(?P<column>\d+)_stats\.txt$''')


def run(h5_filenames, alignment_parameters, alignment_tile_data, all_tile_data, experiment, metadata, make_pdfs):
    if len(h5_filenames) == 0:
        error.fail("There were no HDF5 files to process. "
                   "Either they just don't exist, or you didn't provide the correct path.")
    channel = metadata['alignment_channel']
    experiment_chip = chip.load(metadata['chip_type'])(metadata['ports_on_right'])
    # We use one process per concentration. We could theoretically speed this up since our machine
    # has significantly more cores than the typical number of concentration points, but since it
    # usually finds a result in the first image or two, it's not going to deliver any practical benefits
    num_processes = len(h5_filenames)
    pool = multiprocessing.Pool(num_processes)
    fia = fastqimagealigner.FastqImageAligner(experiment)
    fia.load_reads(alignment_tile_data)

    for h5_filename in h5_filenames:
        base_name = os.path.splitext(h5_filename)[0]
        for directory in (experiment.figure_directory, experiment.results_directory):
            full_directory = os.path.join(directory, base_name)
            if not os.path.exists(full_directory):
                os.makedirs(full_directory)

    with h5py.File(h5_filenames[0]) as first_file:
        grid = GridImages(first_file, channel)
        # find columns/tiles on the left side

        base_column_checker = functools.partial(check_column_for_alignment, channel, alignment_parameters,
                                                metadata['microns_per_pixel'], fia)

        left_end_tiles = dict(get_bounds(pool, h5_filenames, base_column_checker, grid.columns, experiment_chip.left_side_tiles))
        right_end_tiles = dict(get_bounds(pool, h5_filenames, base_column_checker, reversed(grid.columns), experiment_chip.right_side_tiles))

    default_left_tile, default_left_column = decide_default_tiles_and_columns(left_end_tiles)
    default_right_tile, default_right_column = decide_default_tiles_and_columns(right_end_tiles)
    end_tiles = build_end_tiles(h5_filenames, experiment_chip, left_end_tiles, default_left_tile, right_end_tiles,
                                default_right_tile, default_left_column, default_right_column)

    # Iterate over images that are probably inside an Illumina tile, attempt to align them, and if they
    # align, do a precision alignment and write the mapped FastQ reads to disk
    num_processes = multiprocessing.cpu_count()
    log.debug("Aligning all images with %d cores" % num_processes)
    alignment_func = functools.partial(perform_alignment, alignment_parameters, metadata['microns_per_pixel'],
                                       experiment, all_tile_data, make_pdfs, fia)

    pool = multiprocessing.Pool(num_processes)
    pool.map_async(alignment_func, iterate_all_images(h5_filenames, end_tiles, channel), chunksize=96).get(timeout=sys.maxint)
    log.debug("Done aligning!")


def build_end_tiles(h5_filenames, experiment_chip, left_end_tiles, default_left_tile, right_end_tiles,
                    default_right_tile, default_left_column, default_right_column):
    end_tiles = {}
    # Now build up the end tile data structure
    for filename in h5_filenames:
        left_tiles, left_column = left_end_tiles.get(filename, ([default_left_tile], default_left_column))
        right_tiles, right_column = right_end_tiles.get(filename, ([default_right_tile], default_right_column))
        min_column, max_column = min(left_column, right_column), max(left_column, right_column)
        tile_map = experiment_chip.expected_tile_map(left_tiles, right_tiles, min_column, max_column)
        end_tiles[filename] = min_column, max_column, tile_map
    return end_tiles


def run_second_channel(h5_filenames, channel_name, alignment_parameters, all_tile_data, experiment, clargs):
    num_processes = multiprocessing.cpu_count()
    log.debug("Doing second channel alignment of all images with %d cores" % num_processes)
    fastq_image_aligner = fastqimagealigner.FastqImageAligner(experiment)
    fastq_image_aligner.load_reads(all_tile_data)
    second_processor = functools.partial(process_data_image, alignment_parameters, all_tile_data,
                                         clargs.microns_per_pixel, experiment, clargs.make_pdfs,
                                         channel_name, fastq_image_aligner)

    pool = multiprocessing.Pool(num_processes)
    pool.map_async(second_processor,
                   load_aligned_stats_files(h5_filenames, clargs.alignment_channel, experiment)).get(sys.maxint)
    log.debug("Done aligning!")


def extract_rc_info(stats_file):
    match = stats_regex.match(stats_file)
    if match:
        return int(match.group('row')), int(match.group('column'))
    raise ValueError("Invalid stats file: %s" % str(stats_file))


def load_aligned_stats_files(h5_filenames, channel, experiment):
    for h5_filename in h5_filenames:
        base_name = os.path.splitext(h5_filename)[0]
        for f in os.listdir(os.path.join(experiment.results_directory, base_name)):
            if f.endswith('_stats.txt') and channel in f:
                try:
                    row, column = extract_rc_info(f)
                except ValueError:
                    log.warn("Invalid stats file: %s" % str(f))
                    continue
                else:
                    yield h5_filename, base_name, f, row, column


def process_data_image(alignment_parameters, tile_data, um_per_pixel, experiment, make_pdfs, channel, fastq_image_aligner,
                       (h5_filename, base_name, stats_filepath, row, column)):
    with h5py.File(h5_filename) as h5:
        grid = GridImages(h5, channel)
        image = grid.get(row, column)
    sexcat_filepath = os.path.join(base_name, '%s.cat' % image.index)
    stats_filepath = os.path.join(experiment.results_directory, base_name, stats_filepath)
    fastq_image_aligner.set_image_data(image, um_per_pixel)
    fastq_image_aligner.set_sexcat_from_file(sexcat_filepath)
    fastq_image_aligner.alignment_from_alignment_file(stats_filepath)
    try:
        fastq_image_aligner.precision_align_only(min_hits=alignment_parameters.min_hits)
        log.debug("Processed 2nd channel for %s" % image.index)
    except ValueError:
        log.debug("Could not precision align %s" % image.index)
    else:
        write_output(image.index, base_name, fastq_image_aligner, experiment, tile_data, make_pdfs)


def decide_default_tiles_and_columns(end_tiles):
    all_tiles = []
    columns = []
    for filename, (tiles, column) in end_tiles.items():
        for tile in tiles:
            all_tiles.append(tile)
        columns.append(column)
    a, b = Counter(all_tiles).most_common(1), Counter(columns).most_common(1)
    best_tile, best_column = a[0][0], b[0][0]
    return best_tile, best_column


def get_bounds(pool, h5_filenames, base_column_checker, columns, possible_tile_keys):
    end_tiles = Manager().dict()
    for column in columns:
        column_checker = functools.partial(base_column_checker, end_tiles, column, possible_tile_keys)
        pool.map_async(column_checker, h5_filenames).get(sys.maxint)
        if end_tiles:
            return end_tiles
    return False


def check_column_for_alignment(channel, alignment_parameters, um_per_pixel, fia,
                               end_tiles, column, possible_tile_keys, h5_filename):
    base_name = os.path.splitext(h5_filename)[0]
    with h5py.File(h5_filename) as h5:
        grid = GridImages(h5, channel)
        image = grid.get(3, column)
        log.debug("Aligning %s Row 3 Column %d against PhiX" % (base_name, column))
        fia = process_alignment_image(alignment_parameters, base_name, um_per_pixel,
                                      image, possible_tile_keys, deepcopy(fia))
        if fia.hitting_tiles:
            log.debug("%s aligned to at least one tile!" % image.index)
            # because of the way we iterate through the images, if we find one that aligns,
            # we can just stop because that gives us the outermost column of images and the
            # outermost FastQ tile
            end_tiles[h5_filename] = [tile.key for tile in fia.hitting_tiles], image.column


def perform_alignment(alignment_parameters, um_per_pixel, experiment, all_tile_data, make_pdfs, preloaded_fia, image_data):
    # Does a rough alignment, and if that works, does a precision alignment and writes the corrected
    # FastQ reads to disk
    row, column, channel, h5_filename, possible_tile_keys, base_name = image_data
    with h5py.File(h5_filename) as h5:
        grid = GridImages(h5, channel)
        image = grid.get(row, column)
    log.debug("Aligning image from %s. Row: %d, Column: %d " % (base_name, image.row, image.column))
    # first get the correlation to random tiles, so we can distinguish signal from noise
    fia = process_alignment_image(alignment_parameters, base_name, um_per_pixel, image,
                                  possible_tile_keys, deepcopy(preloaded_fia))
    if fia.hitting_tiles:
        # The image data aligned with FastQ reads!
        try:
            fia.precision_align_only(hit_type=('exclusive', 'good_mutual'),
                                     min_hits=alignment_parameters.min_hits)
        except ValueError:
            log.debug("Too few hits to perform precision alignment. Image: %s Row: %d Column: %d " % (base_name, image.row, image.column))
        else:
            write_output(image.index, base_name, fia, experiment, all_tile_data, make_pdfs)
    # The garbage collector takes its sweet time for some reason, so we have to manually delete
    # these objects or memory usage blows up.
    del fia
    del image


def iterate_all_images(h5_filenames, end_tiles, channel):
    # We need an iterator over all images to feed the parallel processes. Since each image is
    # processed independently and in no particular order, we need to return information in addition
    # to the image itself that allow files to be written in the correct place and such
    for h5_filename in h5_filenames:
        base_name = os.path.splitext(h5_filename)[0]
        with h5py.File(h5_filename) as h5:
            grid = GridImages(h5, channel)
            min_column, max_column, tile_map = end_tiles[h5_filename]
            for column in range(min_column, max_column):
                for row in range(grid._height):
                    image = grid.get(row, column)
                    if image is not None:
                        yield row, column, channel, h5_filename, tile_map[image.column], base_name


def load_read_names(file_path):
    # reads a FastQ file with Illumina read names
    with open(file_path) as f:
        tiles = defaultdict(set)
        for line in f:
            lane, tile = line.strip().rsplit(':', 4)[1:3]
            key = 'lane{0}tile{1}'.format(lane, tile)
            tiles[key].add(line.strip())
    del f
    return {key: list(values) for key, values in tiles.items()}


def process_alignment_image(alignment_parameters, base_name, um_per_pixel, image, possible_tile_keys, fia):
    sexcat_fpath = os.path.join(base_name, '%s.cat' % image.index)
    fia.set_image_data(image, um_per_pixel)
    fia.set_sexcat_from_file(sexcat_fpath)
    fia.rough_align(possible_tile_keys,
                    alignment_parameters.rotation_estimate,
                    alignment_parameters.fastq_tile_width_estimate,
                    snr_thresh=alignment_parameters.snr)
    return fia


def write_output(image_index, base_name, fastq_image_aligner, experiment, tile_data, make_pdfs):
    intensity_filepath = os.path.join(experiment.results_directory,
                                      base_name, '{}_intensities.txt'.format(image_index))
    stats_filepath = os.path.join(experiment.results_directory,
                                  base_name, '{}_stats.txt'.format(image_index))
    all_read_rcs_filepath = os.path.join(experiment.results_directory,
                                         base_name, '{}_all_read_rcs.txt'.format(image_index))

    if make_pdfs:
        ax = plotting.plot_all_hits(fastq_image_aligner)
        ax.figure.savefig(os.path.join(experiment.figure_directory, '{}_all_hits.pdf'.format(image_index)))
        ax = plotting.plot_hit_hists(fastq_image_aligner)
        ax.figure.savefig(os.path.join(experiment.figure_directory, '{}_hit_hists.pdf'.format(image_index)))

    fastq_image_aligner.output_intensity_results(intensity_filepath)
    fastq_image_aligner.write_alignment_stats(stats_filepath)
    all_fastq_image_aligner = fastqimagealigner.FastqImageAligner(experiment)
    all_fastq_image_aligner.all_reads_fic_from_aligned_fic(fastq_image_aligner, tile_data)
    all_fastq_image_aligner.write_read_names_rcs(all_read_rcs_filepath)
