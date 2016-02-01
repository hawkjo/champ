import os


class AlignmentParameters(object):
    def __init__(self, base_directory, chip_id, aligned_image_index_offset=0, fq_w_estimate=935.0, min_tile=1, max_tile=19, min_hits=15,
                 objective=60, rotation_estimate=180.0, snr_threshold=1.2, strategy='slow'):
        assert strategy in ('fast', 'slow'), 'Invalid alignment strategy: {strategy}'.format(strategy=strategy)
        self._base_directory = base_directory
        self._chip_id = chip_id
        self._aligned_image_index_offset = aligned_image_index_offset
        self._fq_w_estimate = fq_w_estimate
        self._min_tile = min_tile
        self._max_tile = max_tile
        self._min_hits = min_hits
        self._objective = objective
        self._rotation_estimate = rotation_estimate
        self._snr_threshold = snr_threshold
        self._strategy = strategy

    @property
    def aligning_read_names_filepath(self):
        return self._make_filepath('phiX_mappings/phiX_read_names.txt')

    @property
    def aligned_image_index_offset(self):
        return int(self._aligned_image_index_offset)

    @property
    def all_read_names_filepath(self):
        return self._make_filepath('read_names/all_read_names.txt')

    @property
    def fq_w_est(self):
        return float(self._fq_w_estimate)

    @property
    def max_tile_num(self):
        return 2100 + int(self._max_tile)

    @property
    def min_tile_num(self):
        return 2100 + int(self._min_tile)

    @property
    def min_hits(self):
        return int(self._min_hits)

    @property
    def objective(self):
        return int(self._objective)

    @property
    def chip_id(self):
        return self._chip_id

    @property
    def rotation_estimate(self):
        return float(self._rotation_estimate)

    @property
    def snr_threshold(self):
        return float(self._snr_threshold)

    @property
    def strategy(self):
        return self._strategy

    def _make_filepath(self, filename):
        return '{base_directory}{sep}{chip_id}{sep}{filename}'.format(base_directory=self._base_directory,
                                                                      chip_id=self._chip_id,
                                                                      sep=os.path.sep,
                                                                      filename=filename)
