from chimp.config import AlignmentParameters, Experiment
from chimp.process_nd2_im import process_fig
from chimp.grid import GridImages
import os
import logging
from nd2reader import Nd2


log = logging.getLogger(__name__)


def tile_keys_given_nums(tile_nums):
    return ['lane1tile{0}'.format(tile_num) for tile_num in tile_nums]


def main(clargs):
    nd2_filenames = list(filter(lambda x: x.endswith('.nd2'), os.listdir(os.getcwd())))
    LEFT_TILES = tile_keys_given_nums(range(2101, 2111))
    RIGHT_TILES = tile_keys_given_nums(reversed(range(2111, 2120)))
    experiment = Experiment(clargs.project_name)
    for nd2_filename in nd2_filenames:
        nd2 = Nd2(nd2_filename)
        # CHANNEL OFFSET IS SET TO 1 JUST BECAUSE WE ARE GOING TO REMOVE THIS ENTIRELY
        # WHEN WE SWITCH TO MICROMANAGER
        grid = GridImages(nd2, channel_offset=1)
        log.debug("nd2_filenames %s" % nd2_filenames)
        alignment_parameters = AlignmentParameters(clargs)
        left_tile, left_column = find_end_tile(grid.left_iter(), alignment_parameters, nd2_filename,
                                               LEFT_TILES, experiment)
        right_tile, right_column = find_end_tile(grid.right_iter(), alignment_parameters, nd2_filename,
                                                 RIGHT_TILES, experiment)


def find_end_tile(indexes, alignment_parameters, nd2_filename, possible_tiles, experiment):
    objective = 60
    for index, row, column in indexes:
        nd2 = Nd2(nd2_filename)[index]
        # first get the correlation to random tiles, so we can distinguish signal from noise
        fia = process_fig(alignment_parameters, nd2, nd2_filename, index, objective, possible_tiles, experiment)
        if fia.hitting_tiles:
            print(index, row, column, fia.hitting_tiles)
            return column, fia.hitting_tiles
