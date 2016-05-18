from chimp import constants
import fastqimagealigner
from chimp.grid import GridImages
import functools
import os
import logging
import h5py
from collections import defaultdict
import multiprocessing
import sys

log = logging.getLogger(__name__)


def run(h5_filenames, alignment_parameters, alignment_tile_data, experiment, um_per_pixel, channel):
    end_tiles = {}
    boundary_finder = functools.partial(find_boundary_columns, channel, alignment_parameters,
                                        alignment_tile_data, um_per_pixel, experiment, end_tiles)
    num_processes = 1 # len(h5_filenames)
    pool = multiprocessing.Pool(num_processes)
    log.debug("Finding boundaries with %d processes" % num_processes)
    pool.map_async(boundary_finder, h5_filenames).get(timeout=sys.maxint)

    if not end_tiles:
        log.debug("whoops, no end tiles")
        exit()

    log.debug("Done finding boundaries!")
    print("end_tiles", end_tiles)

    # Iterate over images that are probably inside an Illumina tile, attempt to align them, and if they
    # align, do a precision alignment and write the mapped FastQ reads to disk
    num_processes = multiprocessing.cpu_count()
    log.debug("Aligning all images with %d cores" % num_processes)
    alignment_func = functools.partial(alignment, alignment_parameters, um_per_pixel,
                                       experiment, alignment_tile_data)
    pool = multiprocessing.Pool(num_processes)
    pool.map_async(alignment_func, iterate_all_images(h5_filenames, end_tiles, channel)).get(timeout=sys.maxint)
    log.debug("Done aligning!")


def alignment(alignment_parameters, um_per_pixel, experiment, alignment_tile_data, image_data):
        image, possible_tile_keys, base_name = image_data
        log.debug("Aligning image from %s. Row: %d, Column: %d " % (base_name, image.row, image.column))
        # first get the correlation to random tiles, so we can distinguish signal from noise
        fia = process_fig(alignment_parameters, base_name, alignment_tile_data,  um_per_pixel,
                          experiment, image, possible_tile_keys)
        if fia.hitting_tiles:
            # The image data aligned with FastQ reads!
            fia.precision_align_only(hit_type=('exclusive', 'good_mutual'),
                                     min_hits=alignment_parameters.min_hits)
            write_output(image.index, base_name, fia, experiment, alignment_tile_data)
        # The garbage collector takes its sweet time for some reason, so we have to manually delete
        # these objects or memory usage blow up
        del fia
        del image


def iterate_all_images(h5_filenames, end_tiles, channel):
    for h5_filename in h5_filenames:
        base_name = os.path.splitext(h5_filename)[0]
        with h5py.File(h5_filename) as h5:
            grid = GridImages(h5, channel)
            left_column, right_column, tile_map = end_tiles[h5_filename]
            for image in grid.bounded_iter(left_column, right_column):
                yield image, tile_map[image.column], base_name


def find_boundary_columns(channel, alignment_parameters, alignment_tile_data, um_per_pixel,
                          experiment, end_tiles, h5_filename):
    log.debug("Finding boundary columns")
    # Align image data to FastQ reads and write the aligned FastQ reads to disk
    base_name = os.path.splitext(h5_filename)[0]
    log.debug("fbc base name %s" % base_name)
    with h5py.File(h5_filename) as h5:
        grid = GridImages(h5, channel)

        # We need to call process_fig() several times with almost the same parameters
        figure_processor = functools.partial(process_fig, alignment_parameters, base_name,
                                             alignment_tile_data, um_per_pixel, experiment)

        # Find the outermost columns of image data where we overlap with FastQ tile reads
        # We do this so we can skip any images that are definitely not going to be useful to us
        left_column, right_column, tile_map = find_ends(grid, figure_processor)
        end_tiles[base_name] = left_column, right_column, tile_map


def load_read_names(file_path):
    with open(file_path) as f:
        tiles = defaultdict(set)
        for line in f:
            lane, tile = line.strip().rsplit(':', 4)[1:3]
            key = 'lane{0}tile{1}'.format(lane, tile)
            tiles[key].add(line.strip())
    del f
    return {key: list(values) for key, values in tiles.items()}


def find_ends(grid, figure_processor):
    # Determines which tiles we have image data from, for left and right sides of the chip.
    log.info("Finding end tiles")
    left_side_tiles = range(1, 11)
    right_side_tiles = reversed(range(11, 20))

    right_tiles, right_column = find_end_tile(figure_processor, grid.right_iter(), right_side_tiles)
    left_tiles, left_column = find_end_tile(figure_processor, grid.left_iter(), left_side_tiles)

    # do full alignment for images
    # skip end tile finding for make fast
    tile_map = get_expected_tile_map(left_tiles,
                                     right_tiles,
                                     left_column,
                                     right_column)
    return left_column, right_column, tile_map


def get_expected_tile_map(left_tiles, right_tiles, min_column, max_column):
    # Creates a dictionary that relates each column of microscope images to its expected tile, +/- 1.
    tile_map = defaultdict(list)
    min_tile = min([int(tile[-4:]) for tile in left_tiles])
    max_tile = max([int(tile[-4:]) for tile in right_tiles])
    normalization_factor = float(max_tile - min_tile + 1) / float(max_column - min_column)
    for column in range(min_column, max_column + 1):
        expected_tile_number = min(constants.MISEQ_TILE_COUNT,
                                   max(1, int(round(column * normalization_factor, 0)))) + min_tile - 1
        tile_map[column].append(format_tile_number(expected_tile_number))
        if expected_tile_number > min_tile:
            tile_map[column].append(format_tile_number(expected_tile_number - 1))
        if expected_tile_number < max_tile:
            tile_map[column].append(format_tile_number(expected_tile_number + 1))
    return tile_map


def format_tile_number(number):
    # hardcoding 2000, since we can only physically image side 2 currently.
    return 'lane1tile{0}'.format(2000 + number)


def find_end_tile(figure_processor, images, possible_tiles):
    # Figures out which FastQ tile and column of image data are the furthest to the left or right
    # of the chip. By doing this we don't have to waste time aligning images with tiles that can't
    # possibly go together
    for image in images:
        # first get the correlation to random tiles, so we can distinguish signal from noise
        fia = figure_processor(image, possible_tiles)
        if fia.hitting_tiles:
            log.debug("%s aligned to at least one tile!" % image.index)
            # because of the way we iterate through the images, if we find one that aligns,
            # we can just stop because that gives us the outermost column of images and the
            # outermost FastQ tile
            return fia.hitting_tiles, image.column
        else:
            log.debug("%s did not align to any tiles." % image.index)


def process_fig(alignment_parameters, base_name, tile_data,
                um_per_pixel, experiment, image, possible_tile_keys):
    for directory in (experiment.figure_directory, experiment.results_directory):
        full_directory = os.path.join(directory, base_name)
        if not os.path.exists(full_directory):
            os.makedirs(full_directory)
    sexcat_fpath = os.path.join(base_name, '%s.cat' % image.index)
    fic = fastqimagealigner.FastqImageAligner(experiment)
    fic.load_reads(tile_data)
    fic.set_image_data(image, um_per_pixel)
    fic.set_sexcat_from_file(sexcat_fpath)
    fic.rough_align(possible_tile_keys,
                    alignment_parameters.rotation_estimate,
                    alignment_parameters.fq_w_est,
                    snr_thresh=alignment_parameters.snr_threshold)
    return fic


def process_second_fig(alignment_parameters, base_name, tile_data, um_per_pixel, experiment, image):
    sexcat_filepath = os.path.join(base_name, '%s.cat' % image.index)
    stats_filepath = os.path.join(experiment.results_directory,
                                  base_name,
                                  '{}_stats.txt'.format(image.index))
    fastq_image_aligner = fastqimagealigner.FastqImageAligner(experiment)
    fastq_image_aligner.load_reads(tile_data)
    fastq_image_aligner.set_image_data(image, um_per_pixel)
    fastq_image_aligner.set_sexcat_from_file(sexcat_filepath)
    fastq_image_aligner.alignment_from_alignment_file(stats_filepath)
    fastq_image_aligner.precision_align_only(min_hits=alignment_parameters.min_hits)

    write_output(image.index, base_name, fastq_image_aligner, experiment, tile_data)


def write_output(image_index, base_name, fastq_image_aligner, experiment, tile_data):
    intensity_filepath = os.path.join(experiment.results_directory,
                                      base_name, '{}_intensities.txt'.format(image_index))
    stats_filepath = os.path.join(experiment.results_directory,
                                  base_name, '{}_stats.txt'.format(image_index))
    all_read_rcs_filepath = os.path.join(experiment.results_directory,
                                         base_name, '{}_all_read_rcs.txt'.format(image_index))

    fastq_image_aligner.output_intensity_results(intensity_filepath)
    fastq_image_aligner.write_alignment_stats(stats_filepath)
    all_fastq_image_aligner = fastqimagealigner.FastqImageAligner(experiment)
    all_fastq_image_aligner.all_reads_fic_from_aligned_fic(fastq_image_aligner, tile_data)
    all_fastq_image_aligner.write_read_names_rcs(all_read_rcs_filepath)
