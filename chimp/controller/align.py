from chimp.config import AlignmentParameters, Experiment
from chimp.process_nd2_im import process_fig, precision_process_fig, write_output
from chimp.grid import GridImages
import os
import logging
from nd2reader import Nd2
from collections import defaultdict
from chimp.model import constants


log = logging.getLogger(__name__)


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


def main(clargs):
    nd2_filenames = list(filter(lambda x: x.endswith('.nd2'), os.listdir(os.getcwd())))
    LEFT_TILES = tile_keys_given_nums(range(2101, 2111))
    RIGHT_TILES = tile_keys_given_nums(reversed(range(2111, 2120)))
    experiment = Experiment(clargs.project_name)
    objective = 60
    for nd2_filename in nd2_filenames:
        nd2 = Nd2(nd2_filename)
        # CHANNEL OFFSET IS SET TO 1 JUST BECAUSE WE ARE GOING TO REMOVE THIS ENTIRELY
        # WHEN WE SWITCH TO MICROMANAGER
        grid = GridImages(nd2, channel_offset=1)
        alignment_parameters = AlignmentParameters(clargs)
        left_tiles, left_column = find_end_tile(grid.left_iter(), alignment_parameters, nd2_filename,
                                                LEFT_TILES, experiment, objective)
        left_tile = min([int(tile.key[-4:]) for tile in left_tiles])
        right_tiles, right_column = find_end_tile(grid.right_iter(), alignment_parameters, nd2_filename,
                                                  RIGHT_TILES, experiment, objective)
        right_tile = min([int(tile.key[-4:]) for tile in right_tiles])
        # do full alignment for images

        tile_map = get_expected_tile_map(left_tile - 2100, right_tile - 2100, left_column, right_column)
        for index, row, column in grid.bounded_iter(left_column, right_column):
            log.debug("Aligning image #%d from %s" % (index, nd2_filename))
            image = nd2[index]
            print("tile map", tile_map[column])
            tile_numbers = (2100 + tile for tile in tile_map[column])
            possible_tiles = tile_keys_given_nums(tile_numbers)
            print("possible tiles", possible_tiles)
            # first get the correlation to random tiles, so we can distinguish signal from noise
            fia = process_fig(alignment_parameters,
                              image,
                              nd2_filename,
                              index,
                              objective,
                              possible_tiles,
                              experiment)
            if fia.hitting_tiles:
                print("About to precision align")
                precision_process_fig(fia, alignment_parameters)
                write_output(index, fia, experiment, alignment_parameters)
            else:
                print("#%d did not align to %s" % (index, " ".join(possible_tiles)))


def find_end_tile(indexes, alignment_parameters, nd2_filename, possible_tiles, experiment, objective):
    nd2 = Nd2(nd2_filename)
    for index, row, column in indexes:
        image = nd2[index]
        # first get the correlation to random tiles, so we can distinguish signal from noise
        fia = process_fig(alignment_parameters,
                          image,
                          nd2_filename,
                          index,
                          objective,
                          possible_tiles,
                          experiment)
        if fia.hitting_tiles:
            return fia.hitting_tiles, column
