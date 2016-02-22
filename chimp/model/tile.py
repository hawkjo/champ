import numpy as np
from model.constants import FASTQ_TILE_WIDTH


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

    def _load(self, um_per_pixel, read_data):
        tile_shapes = []

        for region, fastq_reads in read_data.items():
            rcs = np.array([(r.row, r.column) for r in fastq_reads])
            tile_shapes.append(rcs.max(axis=0) + 1)
            self._tiles[region] = Tile(region, rcs)

        all_data = np.concatenate([tile.rcs for tile in self._tiles.values()])
        x_min, y_min = all_data.min(axis=0)
        x_max, y_max = all_data.max(axis=0)
        self._offset = np.array([-x_min, -y_min])
        self._scale = (FASTQ_TILE_WIDTH / (x_max - x_min)) / um_per_pixel
        self._largest_tile_shape = np.array(tile_shapes).max(axis=0)

    def fft(self, tile_number, lane=1, side=2):
        key = (lane, '%d1%2d' % (side, tile_number))
        fft = self._ffts.get(key)
        if fft is None:
            tile = self._tiles[key]
            self._ffts[key] = tile.fft(self._scale, self._offset, self._largest_tile_shape)
        return self._ffts[key]


class Tile(object):
    """ Wraps fastq tile coordinates """
    def __init__(self, region, rcs):
        # not sure if tile needs to know its region anymore
        self.region = region
        self._rcs = rcs

    @property
    def rcs(self):
        return self._rcs

    def fft(self, scale, offset, largest_tile_shape):
        tile_image = self._flipped()._image(scale, offset, largest_tile_shape)
        return np.fft.fft2(tile_image)

    def _flipped(self):
        flip_180 = np.array([[-1, 0], [0, -1]])
        flipped_rcs = np.dot(self.rcs, flip_180)
        flipped_rcs -= np.tile(flipped_rcs.min(axis=0), (flipped_rcs.shape[0], 1))
        return Tile(self.region, flipped_rcs)

    def _image(self, scale, offset, largest_tile_shape):
        new_rcs = scale * (self.rcs + np.tile(offset, (self.rcs.shape[0], 1)))
        image = np.zeros(largest_tile_shape)
        image[new_rcs.astype(np.int)[:, 0], new_rcs.astype(np.int)[:, 1]] = 1
        return image
