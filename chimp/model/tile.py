import numpy as np
from model.constants import FASTQ_TILE_WIDTH


class TileManager(object):
    """ Provides access to Tiles. """
    def __init__(self):
        self._tiles = {}

    def load_for_alignment(self, um_per_pixel, read_data):
        tile_shapes = []

        for region, fastq_reads in read_data.items():
            rcs = np.array([(r.row, r.column) for r in fastq_reads])
            flip_180 = np.array([[-1, 0], [0, -1]])
            flipped_rcs = np.dot(rcs, flip_180)
            flipped_rcs -= np.tile(flipped_rcs.min(axis=0), (flipped_rcs.shape[0], 1))
            tile_shapes.append(flipped_rcs.max(axis=0) + 1)
            self._tiles[region] = Tile(region, flipped_rcs)

        all_data = np.concatenate([tile.rcs for tile in self._tiles.values()])
        x_min, y_min = all_data.min(axis=0)
        x_max, y_max = all_data.max(axis=0)
        offset = np.array([-x_min, -y_min])
        scale = (FASTQ_TILE_WIDTH / (x_max - x_min)) / um_per_pixel
        scaled_maxes = scale * np.array([x_max - x_min, y_max - y_min])
        scaled_dims = np.array(tile_shapes).max(axis=0)
        # JIM: This is half implemented. you're trying to account for all the stuff that gets set in
        # set_fastq_image_data() and rotate_data() method or whatever.
        # Once you figure out all those base transformations, you'll
        # apply them to the tile data BEFORE giving it to the tile.


class Tile(object):
    """ Wraps fastq tile coordinates """
    def __init__(self, region, rcs):
        # not sure if tile needs to know its region anymore
        self.region = region
        self._rcs = rcs

    @property
    def rcs(self):
        return self._rcs
