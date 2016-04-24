from collections import defaultdict
from chimp.model import constants
from chimp.process_nd2_im import process_fig, write_output
from chimp.grid import GridImages
import functools
import time
import os
import logging
import h5py


log = logging.getLogger(__name__)


def run(alignment_parameters, alignment_tile_data, all_tile_data, experiment, objective, h5_filename):
    # Align image data to FastQ reads and write the aligned FastQ reads to disk
    print("align_image_data!!! %s" % h5_filename)
    base_name = os.path.splitext(h5_filename)[0]
    h5 = h5py.File(h5_filename)
    grid = GridImages(h5, alignment_parameters.alignment_channel)

    # We need to call process_fig() several times with almost the same parameters
    figure_processor = functools.partial(process_fig, alignment_parameters, base_name,
                                         alignment_tile_data, objective, experiment)

    # Find the outermost columns of image data where we overlap with FastQ tile reads
    # We do this so we can skip any images that are definitely not going to be useful to us
    left_column, right_column, tile_map = find_ends(grid, figure_processor)

    # Iterate over images that are probably inside an Illumina tile, attempt to align them, and if they
    # align, do a precision alignment and write the mapped FastQ reads to disk
    for image in grid.bounded_iter(left_column, right_column):
        start = time.time()
        log.debug("Aligning image from %s. Row: %d, Column: %d " % (h5_filename, image.row, image.column))
        # first get the correlation to random tiles, so we can distinguish signal from noise
        fia = figure_processor(image, tile_map[image.column])
        if fia.hitting_tiles:
            # The image data aligned with FastQ reads!
            fia.precision_align_only(hit_type=('exclusive', 'good_mutual'),
                                     min_hits=alignment_parameters.min_hits)
            write_output(image.index, base_name, fia, experiment, all_tile_data)
        print("%s, row %s column %s took %s seconds to align" % (h5_filename, image.row, image.column, (time.time() - start)))
        # The garbage collector takes its sweet time for some reason, so we have to manually delete
        # these objects or memory usage blow up
        del fia
        del image


def find_ends(grid, figure_processor):
    # Determines which tiles we have image data from, for left and right sides of the chip.
    log.info("Finding end tiles")
    left_side_tiles = range(1, 11)
    right_side_tiles = reversed(range(11, 20))

    left_tiles, left_column = find_end_tile(figure_processor, grid.left_iter(), left_side_tiles)
    right_tiles, right_column = find_end_tile(figure_processor, grid.right_iter(), right_side_tiles)

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
            # because of the way we iterate through the images, if we find one that aligns,
            # we can just stop because that gives us the outermost column of images and the
            # outermost FastQ tile
            return fia.hitting_tiles, image.column
