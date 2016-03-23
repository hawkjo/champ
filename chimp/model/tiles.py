import numpy as np
from chimp.model.constants import FASTQ_TILE_WIDTH
import logging
from chimp import padding
from collections import deque


log = logging.getLogger(__name__)


class TileManager(object):
    def __init__(self, tile_data, scale, offset, image_height, image_width, cache_size):
        self._tiles = {}
        self._fft_conjugate = {}
        self._padding_shape = image_height, image_width
        self._conjugate_cache = {}
        self._cache_tile_numbers = deque()
        self._cache_size = cache_size
        for region, rcs in tile_data.items():
            self._tiles[region] = Tile(flip_coordinates(rcs), scale, offset)

    def _calculate_fft_conjugate(self, image, image_height, image_width):
        """
        Precompute the conjugate of the FFT of the tile points since it's slow and we need to reuse this result many times.

        """
        padding_dimension = padding.calculate_pad_size(image.shape[0], image.shape[1], image_height, image_width)
        padded_image = padding.pad_image(image, padding_dimension)
        tile_fft = np.fft.fft2(padded_image)
        return np.conj(tile_fft)

    def get(self, tile_number, lane=1, side=2):
        return self._tiles[(lane, side, tile_number)]

    def fft_conjugate(self, tile_number, lane=1, side=2):
        region = lane, side, tile_number
        if region not in self._conjugate_cache:
            if len(self._conjugate_cache) >= self._cache_size:
                # cache is full, so evict the oldest member
                oldest = self._cache_tile_numbers.pop()
                del self._conjugate_cache[(lane, side, oldest)]
            self._conjugate_cache[region] = self._calculate_fft_conjugate(self._tiles[region].image,
                                                                          self._padding_shape[0],
                                                                          self._padding_shape[1])
            self._cache_tile_numbers.appendleft(tile_number)
        else:
            # deprioritize this tile number for eviction
            self._cache_tile_numbers.remove(tile_number)
            self._cache_tile_numbers.appendleft(tile_number)
        return self._conjugate_cache[region]


class Tile(object):
    """ Wraps fastq tile coordinates """
    def __init__(self, rcs, scale, offset):
        self._rcs = rcs
        self._scale = scale
        self._offset = offset

    @property
    def rcs(self):
        return self._rcs

    @property
    def normalized_rcs(self):
        return self._scale * (self.rcs + np.tile(self._offset, (self.rcs.shape[0], 1)))

    @property
    def image(self):
        new_rcs = self.normalized_rcs
        image = np.zeros(new_rcs.max(axis=0) + 1)
        image[new_rcs.astype(np.int)[:, 0], new_rcs.astype(np.int)[:, 1]] = 1
        return image


def flip_coordinates(rcs):
    flip_180 = np.array([[-1, 0], [0, -1]])
    flipped_rcs = np.dot(rcs, flip_180)
    flipped_rcs -= np.tile(flipped_rcs.min(axis=0), (flipped_rcs.shape[0], 1))
    return flipped_rcs


def load_tile_manager(um_per_pixel, image_height, image_width, read_data, cache_size):
    tile_shapes = []
    tile_data = {}
    for region, fastq_reads in read_data.items():
        # fastq_reads is a list of FastqRead
        rcs = np.array([(r.column, r.row) for r in fastq_reads])
        tile_shapes.append(rcs.max(axis=0) + 1)
        tile_data[region] = rcs

    all_data = np.concatenate([rcs for rcs in tile_data.values()])
    x_min, y_min = all_data.min(axis=0)
    x_max, y_max = all_data.max(axis=0)
    scale = (FASTQ_TILE_WIDTH / (x_max - x_min)) / um_per_pixel
    offset = np.array([-x_min, -y_min])
    return TileManager(tile_data, scale, offset, image_height, image_width, cache_size)
