import numpy as np
from chimp.model.constants import FASTQ_TILE_WIDTH
import logging

log = logging.getLogger(__name__)


class TileManager(object):
    def __init__(self, tile_data, scale, offset):
        self._tiles = {}
        for region, rcs in tile_data.items():
            self._tiles[region] = Tile(rcs, scale, offset)

    def get(self, tile_number, lane=1, side=2):
        key = lane, side, tile_number
        return self._tiles[key]


class Tile(object):
    """ Wraps fastq tile coordinates """
    def __init__(self, rcs, scale, offset):
        # not sure if tile needs to know its region anymore
        self._rcs = rcs
        self._scale = scale
        self._offset = offset

    @property
    def rcs(self):
        return self._rcs

    @property
    def image(self):
        new_rcs = self._scale * (self.rcs + np.tile(self._offset, (self.rcs.shape[0], 1)))
        image = np.zeros(new_rcs.max(axis=0) + 1)
        image[new_rcs.astype(np.int)[:, 0], new_rcs.astype(np.int)[:, 1]] = 1
        return image


def flip_coordinates(rcs):
    flip_180 = np.array([[-1, 0], [0, -1]])
    flipped_rcs = np.dot(rcs, flip_180)
    flipped_rcs -= np.tile(flipped_rcs.min(axis=0), (flipped_rcs.shape[0], 1))
    return flipped_rcs


def load_tile_manager(um_per_pixel, read_data):
    tile_shapes = []
    tile_data = {}
    for region, fastq_reads in read_data.items():
        # fastq_reads is a list of FastqRead
        rcs = np.array([(r.row, r.column) for r in fastq_reads])
        tile_shapes.append(rcs.max(axis=0) + 1)
        tile_data[region] = flip_coordinates(rcs)

    all_data = np.concatenate([rcs for rcs in tile_data.values()])
    print("all data len", len(all_data))
    x_min, y_min = all_data.min(axis=0)
    x_max, y_max = all_data.max(axis=0)
    print("x_min, y_min", x_min, y_min)
    print("x_max, y_max", x_max, y_max)
    scale = (FASTQ_TILE_WIDTH / (x_max - x_min)) / um_per_pixel
    offset = np.array([-x_min, -y_min])
    print("scale", scale)
    print("offset", offset)
    exit()
    return TileManager(tile_data, scale, offset)
