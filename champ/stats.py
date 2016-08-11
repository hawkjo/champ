import numpy as np
import yaml


class AlignmentStats(object):
    """
    The transformations needed to align an image to FASTQ data, and some data to measure the quality of the alignment

    """
    def __init__(self):
        self._data = {}

    def from_file(self, fh):
        self._data = yaml.load(fh)
        return self

    def from_data(self, tile_keys, scalings, tile_widths, rotations, rc_offsets, hits):
        assert len(tile_keys) == len(scalings) == len(tile_widths) == len(rotations) == len(rc_offsets)
        self._data['tile_keys'] = tile_keys
        self._data['scalings'] = scalings
        self._data['tile_widths'] = tile_widths
        self._data['rotations'] = [rotation * np.pi / 180 for rotation in rotations]
        self._data['rc_offsets'] = rc_offsets
        self._data['hits'] = hits
        return self

    @property
    def score(self):
        # A somewhat arbitrary metric to determine if one alignment is better than another
        score = self._data['hits']['exclusive'] + self._data['hits']['good_mutual']
        print("SCORE", score)
        return score

    def __iter__(self):
        for tile_key, scaling, tile_width, rotation, rc_offset in zip(self._data['tile_keys'],
                                                                      self._data['scalings'],
                                                                      self._data['tile_widths'],
                                                                      self._data['rotations'],
                                                                      self._data['rc_offsets']):
            # The number of hits is an aggregation of all tiles that align to the image, which is why it doesn't
            # change between each iteration.
            yield tile_key, scaling, tile_width, rotation, rc_offset, self._data['hits']

    def __repr__(self):
        return yaml.dump(self._data)
