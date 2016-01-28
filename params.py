class AlignmentParameters(object):
    def __init__(self, lines):
        self._data = {}
        for line in map(str.strip, lines):
            if not line:
                continue
            name, value = line.split()
            self._data[name] = value

        self._snr_threshold = self._data.get('snr_thresh', 1.2)
        self._min_hits = self._data.get('min_hits', 15)
        self._strategy = self._data.get('strategy', 'slow')
        assert self._strategy in ['fast', 'slow'], 'Invalid alignment strategy: {strategy}'.format(strategy=self._strategy)

    @property
    def max_tile_num(self):
        return int(self._data['max_tile_num'])

    @property
    def strategy(self):
        return self._strategy

    @property
    def min_tile_num(self):
        return int(self._data['min_tile_num'])

    @property
    def project_name(self):
        return self._data['project_name']

    @property
    def aligning_read_names_fpath(self):
        return self._data['aligning_read_names_fpath']

    @property
    def all_read_names_fpath(self):
        return self._data['all_read_names_fpath']

    @property
    def objective(self):
        return int(self._data['objective'])

    @property
    def aligned_im_idx_offset(self):
        return int(self._data['aligned_im_idx_offset'])

    @property
    def snr_threshold(self):
        return float(self._snr_threshold)

    @property
    def rotation_estimate(self):
        return float(self._data['rotation_est'])

    @property
    def fq_w_est(self):
        return float(self._data['fq_w_est'])

    @property
    def min_hits(self):
        return int(self._min_hits)


def get_align_params(align_param_fpath):
    with open(align_param_fpath) as lines:
        return AlignmentParameters(lines)
