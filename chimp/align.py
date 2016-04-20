from collections import defaultdict
from chimp.model import constants


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
