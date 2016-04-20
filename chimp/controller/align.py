from chimp.config import AlignmentParameters, Experiment
from chimp.process_nd2_im import process_fig, write_output
from chimp.grid import GridImages
from chimp import reads
import os
import logging
from nd2reader import Nd2
import time
import sys
from chimp import align
import multiprocessing
import functools

log = logging.getLogger(__name__)


def main(clargs):
    nd2_filenames = list(filter(lambda x: x.endswith('.nd2'), os.listdir(os.getcwd())))
    experiment = Experiment(clargs.project_name)
    objective = 60
    alignment_parameters = AlignmentParameters(clargs)
    all_tile_data = reads.get_read_names(alignment_parameters.all_read_names_filepath)
    phix_tile_data = reads.get_read_names(os.path.join(experiment.project_name,
                                                       alignment_parameters.aligning_read_names_filepath))
    processes = min(len(nd2_filenames), multiprocessing.cpu_count())
    log.info("Using %d processes for alignment" % processes)
    pool = multiprocessing.Pool(processes=processes)
    alignment_function = functools.partial(align_image_data,
                                           alignment_parameters,
                                           phix_tile_data,
                                           all_tile_data,
                                           experiment,
                                           objective)
    log.info("Starting alignment")
    pool.map_async(alignment_function, nd2_filenames).get(sys.maxint)
    log.info("Alignment complete")


def align_image_data(alignment_parameters, alignment_tile_data, all_tile_data, experiment, objective, nd2_filename):
    log.info("align_image_data!!!")
    return True
    log.info("I am a blasphemy unto the Lord")
    LEFT_TILES = align.tile_keys_given_nums(range(2101, 2111))
    RIGHT_TILES = align.tile_keys_given_nums(reversed(range(2111, 2120)))
    base_name = os.path.splitext(nd2_filename)[0]
    nd2 = Nd2(nd2_filename)
    # CHANNEL OFFSET IS SET TO 1 JUST BECAUSE WE ARE GOING TO REMOVE THIS ENTIRELY
    # WHEN WE SWITCH TO MICROMANAGER
    grid = GridImages(nd2, channel_offset=1)

    log.info("Finding end tiles")
    left_tiles, left_column = align.find_end_tile(grid.left_iter(),
                                                  alignment_parameters,
                                                  base_name,
                                                  alignment_tile_data,
                                                  LEFT_TILES,
                                                  experiment,
                                                  objective)
    left_tile = min([int(tile.key[-4:]) for tile in left_tiles])
    right_tiles, right_column = align.find_end_tile(grid.right_iter(),
                                                    alignment_parameters,
                                                    base_name,
                                                    alignment_tile_data,
                                                    RIGHT_TILES,
                                                    experiment,
                                                    objective)
    right_tile = min([int(tile.key[-4:]) for tile in right_tiles])
    # do full alignment for images
    # skip end tile finding for make fast
    tile_map = align.get_expected_tile_map(left_tile - 2100,
                                           right_tile - 2100,
                                           left_column,
                                           right_column)
    for index, row, column in grid.bounded_iter(left_column, right_column):
        start = time.time()
        log.debug("Aligning image (%d, %d) from %s" % (row, column, nd2_filename))
        image = nd2[index]
        tile_numbers = (2100 + tile for tile in tile_map[column])
        possible_tiles = align.tile_keys_given_nums(tile_numbers)
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
        print(time.time() - start)
        del fia
        del image
