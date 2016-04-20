from collections import defaultdict
from chimp.model import constants
from chimp.process_nd2_im import process_fig, write_output
from chimp.grid import GridImages
from nd2reader import Nd2
import time
import os
import logging

log = logging.getLogger(__name__)


def run(alignment_parameters, alignment_tile_data, all_tile_data, experiment, objective, nd2_filename):
    print("align_image_data!!! %s" % nd2_filename)
    LEFT_TILES = tile_keys_given_nums(range(2101, 2111))
    RIGHT_TILES = tile_keys_given_nums(reversed(range(2111, 2120)))
    base_name = os.path.splitext(nd2_filename)[0]
    nd2 = Nd2(nd2_filename)
    # CHANNEL OFFSET IS SET TO 1 JUST BECAUSE WE ARE GOING TO REMOVE THIS ENTIRELY
    # WHEN WE SWITCH TO MICROMANAGER
    grid = GridImages(nd2, channel_offset=1)
    log.info("Finding end tiles")
    left_tiles, left_column = find_end_tile(grid.left_iter(),
                                            alignment_parameters,
                                            base_name,
                                            alignment_tile_data,
                                            LEFT_TILES,
                                            experiment,
                                            objective)
    left_tile = min([int(tile.key[-4:]) for tile in left_tiles])
    right_tiles, right_column = find_end_tile(grid.right_iter(),
                                              alignment_parameters,
                                              base_name,
                                              alignment_tile_data,
                                              RIGHT_TILES,
                                              experiment,
                                              objective)
    right_tile = min([int(tile.key[-4:]) for tile in right_tiles])
    # do full alignment for images
    # skip end tile finding for make fast
    tile_map = get_expected_tile_map(left_tile - 2100,
                                     right_tile - 2100,
                                     left_column,
                                     right_column)
    for index, row, column in grid.bounded_iter(left_column, right_column):
        start = time.time()
        log.debug("Aligning image from %s. Row: %d, Column: %d " % (nd2_filename, row, column))
        image = nd2[index]
        tile_numbers = (2100 + tile for tile in tile_map[column])
        possible_tiles = tile_keys_given_nums(tile_numbers)
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
            fia.precision_align_only(hit_type=('exclusive', 'good_mutual'),
                                     min_hits=alignment_parameters.min_hits)
            write_output(index, base_name, fia, experiment, all_tile_data)
        print("%s, row %s column %s took %s seconds to align" % (nd2_filename, row, column, (time.time() - start)))
        del fia
        del image


def tile_keys_given_nums(tile_nums):
    return ['lane1tile{0}'.format(tile_num) for tile_num in tile_nums]


def get_expected_tile_map(min_tile, max_tile, min_column, max_column):
    """
    Creates a dictionary that relates each column of microscope images to its expected tile, +/- 1.

    """
    tile_map = defaultdict(list)
    normalization_factor = float(max_tile - min_tile + 1) / float(max_column - min_column)
    for column in range(min_column, max_column + 1):
        expected_tile = min(constants.MISEQ_TILE_COUNT,
                            max(1, int(round(column * normalization_factor, 0)))) + min_tile - 1
        tile_map[column].append(expected_tile)
        if expected_tile > min_tile:
            tile_map[column].append(expected_tile - 1)
        if expected_tile < max_tile:
            tile_map[column].append(expected_tile + 1)
    return tile_map


def find_end_tile(indexes, alignment_parameters, base_name, alignment_tile_data, possible_tiles, experiment, objective):
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
