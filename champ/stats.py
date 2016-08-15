import numpy as np
import yaml


class AlignmentStats(object):
    """
    The transformations needed to align an image to FASTQ data, and some data to measure the quality of the alignment

    """
    def __init__(self):
        self._data = {}

    def from_file(self, fh):
        assert not isinstance(fh, str)
        self._data = yaml.load(fh)
        self._validate_data()
        return self

    def _validate_data(self):
        if not len(self._data['tile_keys']) == len(self._data['scalings']) == len(self._data['tile_widths']) == len(self._data['rotations']) == len(self._data['rc_offsets']):
            raise ValueError("Corrupt or invalid AlignmentStats file")

    def from_data(self, tile_keys, scalings, tile_widths, rotations, rc_offsets, hits):
        self._data['tile_keys'] = tile_keys
        self._data['scalings'] = scalings
        self._data['tile_widths'] = tile_widths
        self._data['rotations'] = [rotation * np.pi / 180 for rotation in rotations]
        self._data['rc_offsets'] = rc_offsets
        self._data['hits'] = hits
        self._validate_data()
        return self

    @property
    def score(self):
        # A somewhat arbitrary metric to determine if one alignment is better than another
        return self._data['hits']['exclusive'] + self._data['hits']['good_mutual']

    def __iter__(self):
        for tile_key, scaling, tile_width, rotation, rc_offset in zip(self._data['tile_keys'],
                                                                      self._data['scalings'],
                                                                      self._data['tile_widths'],
                                                                      self._data['rotations'],
                                                                      self._data['rc_offsets']):
            # The number of hits is an aggregation of all tiles that align to the image, which is why it doesn't
            # change between each iteration.
            yield tile_key, scaling, tile_width, rotation, rc_offset, self._data['hits']

    @property
    def serialized(self):
        return yaml.dump(dict(self._data))
