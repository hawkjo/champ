import numpy as np
from model.constants import FASTQ_TILE_WIDTH
import logging

log = logging.getLogger(__name__)


class TileManager(object):
    """ Provides access to Tiles. """
    def __init__(self, um_per_pixel, read_data):
        self._tiles = {}
        self._um_per_pixel = um_per_pixel
        self._largest_tile_shape = (0, 0)
        self._offset = np.array([0, 0])
        self._scale = 1.0
        self._load(um_per_pixel, read_data)
        self._ffts = {}

    def tile(self, tile_number, lane=1, side=2):
        key = lane, side, tile_number
        return self._tiles[key]

    def fft(self, tile_number, lane=1, side=2):
        key = lane, side, tile_number
        fft = self._ffts.get(key)
        if fft is None:
            # cache results
            tile = self._tiles[key]
            self._ffts[key] = tile.fft(self._scale, self._offset, self._largest_tile_shape)
        return self._ffts[key]

    def _flip_coordinates(self, rcs):
        flip_180 = np.array([[-1, 0], [0, -1]])
        flipped_rcs = np.dot(rcs, flip_180)
        flipped_rcs -= np.tile(flipped_rcs.min(axis=0), (flipped_rcs.shape[0], 1))
        return flipped_rcs

    def _load(self, um_per_pixel, read_data):
        tile_shapes = []
        tile_data = {}
        for region, fastq_reads in read_data.items():
            # fastq is a list of FastqRead
            rcs = np.array([(r.row, r.column) for r in fastq_reads])
            tile_shapes.append(rcs.max(axis=0) + 1)
            tile_data[region] = rcs

        all_data = np.concatenate([rcs for rcs in tile_data.values()])
        x_min, y_min = all_data.min(axis=0)
        x_max, y_max = all_data.max(axis=0)
        offset = np.array([-x_min, -y_min])
        log.debug("Offset: %s" % offset)
        scale = (FASTQ_TILE_WIDTH / (x_max - x_min)) / um_per_pixel
        log.debug("Scale: %s" % scale)
        largest_tile_shape = np.array(tile_shapes).max(axis=0)
        log.debug("Largest tile shape: %s" % largest_tile_shape)
        for region, rcs in tile_data.items():
            self._tiles[region] = Tile(region, self._flip_coordinates(rcs), scale, offset, largest_tile_shape)


class Tile(object):
    """ Wraps fastq tile coordinates """
    def __init__(self, region, rcs, scale, offset, tile_shape):
        # not sure if tile needs to know its region anymore
        self.region = region
        self._rcs = rcs
        self._scale = scale
        self._offset = offset
        self._tile_shape = tile_shape

    @property
    def rcs(self):
        return self._rcs

    @property
    def fft(self):
        return np.fft.fft2(self.image)

    @property
    def image(self):
        new_rcs = self._scale * (self.rcs + np.tile(self._offset, (self.rcs.shape[0], 1)))
        image = np.zeros(self._tile_shape)
        image[new_rcs.astype(np.int)[:, 0], new_rcs.astype(np.int)[:, 1]] = 1
        return image
