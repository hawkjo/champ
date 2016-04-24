from collections import defaultdict
from chimp.model import constants
from chimp.process_nd2_im import process_fig, write_output
from chimp.grid import GridImages
from nd2reader import Nd2
import functools
import time
import os
import logging

log = logging.getLogger(__name__)


def run(alignment_parameters, alignment_tile_data, all_tile_data, experiment, objective, nd2_filename):
    print("align_image_data!!! %s" % nd2_filename)
    base_name = os.path.splitext(nd2_filename)[0]
    nd2 = Nd2(nd2_filename)
    # CHANNEL OFFSET IS SET TO 1 JUST BECAUSE WE ARE GOING TO REMOVE THIS ENTIRELY
    # WHEN WE SWITCH TO MICROMANAGER
    grid = GridImages(nd2, channel_offset=1)
    end_tile_finder = functools.partial(find_end_tile, alignment_parameters, base_name,
                                        alignment_tile_data, experiment, objective)
    left_column, right_column, tile_map = find_ends(grid, end_tile_finder)

    for index, row, column in grid.bounded_iter(left_column, right_column):
        start = time.time()
        log.debug("Aligning image from %s. Row: %d, Column: %d " % (nd2_filename, row, column))
        image = nd2[index]
        # first get the correlation to random tiles, so we can distinguish signal from noise
        fia = process_fig(alignment_parameters,
                          image,
                          base_name,
                          alignment_tile_data,
                          index,
                          objective,
                          tile_map[column],
                          experiment)
        if fia.hitting_tiles:
            fia.precision_align_only(hit_type=('exclusive', 'good_mutual'),
                                     min_hits=alignment_parameters.min_hits)
            write_output(index, base_name, fia, experiment, all_tile_data)
        print("%s, row %s column %s took %s seconds to align" % (nd2_filename, row, column, (time.time() - start)))
        del fia
        del image


def find_ends(grid, end_tile_finder):
    """ Determines which tiles we have image data from, for left and right sides of the chip. """
    log.info("Finding end tiles")
    left_tiles = range(1, 11)
    right_tiles = reversed(range(11, 20))

    left_tiles, left_column = end_tile_finder(grid.left_iter(), left_tiles)
    right_tiles, right_column = end_tile_finder(grid.right_iter(), right_tiles)

    # do full alignment for images
    # skip end tile finding for make fast
    tile_map = get_expected_tile_map(left_tiles,
                                     right_tiles,
                                     left_column,
                                     right_column)
    return left_column, right_column, tile_map


def get_expected_tile_map(left_tiles, right_tiles, min_column, max_column):
    """
    Creates a dictionary that relates each column of microscope images to its expected tile, +/- 1.

    """
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
    return 'lane1tile{0}'.format(2000 + number)


def find_end_tile(alignment_parameters, base_name, alignment_tile_data, experiment, objective, indexes, possible_tiles):
    nd2 = Nd2(base_name + ".nd2")
    for index, row, column in indexes:
        image = nd2[index]
        # first get the correlation to random tiles, so we can distinguish signal from noise
        fia = process_fig(alignment_parameters,
                          image,
                          base_name,
                          alignment_tile_data,
                          index,
                          objective,
                          possible_tiles,
                          experiment)
        if fia.hitting_tiles:
            return fia.hitting_tiles, column

