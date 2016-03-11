import os


class Experiment(object):
    def __init__(self, base_dir):
        self._base_dir = base_dir

    @property
    def figure_directory(self):
        return os.path.join(self._base_dir, 'figs')

    @property
    def data_directory(self):
        return os.path.join(self._base_dir, 'data')

    @property
    def fourier_data_directory(self):
        return os.path.join(self.data_directory, 'from_fourierseq')

    @property
    def results_directory(self):
        return os.path.join(self._base_dir, 'results')

    def get_sexcat_path(self, nd2_name, image_index):
        return os.path.join(self._base_dir, nd2_name, '%d.cat' % image_index)


class AlignmentParameters(object):
    """ Parses user-provided alignment parameters and provides a default in case no value was given. """
    def __init__(self, command_line_args, base_directory, presumptive_chip_id):
        self._args = command_line_args
        self._base_directory = base_directory
        self._presumptive_chip_id = presumptive_chip_id

    @property
    def aligning_read_names_filepath(self):
        return self._make_filepath('phiX_mappings/phiX_read_names.txt')

    @property
    def aligned_image_index_offset(self):
        return int(self._args.get('--aligned_image_index_offset') or 0)

    @property
    def all_read_names_filepath(self):
        return self._make_filepath('read_names/all_read_names.txt')

    @property
    def base_directory(self):
        return self._base_directory

    @property
    def fq_w_est(self):
        return float(self._args.get('--fq_w_estimate') or 935.0)

    @property
    def max_tile_num(self):
        return 2100 + int(self._args.get('--max_tile') or 19)

    @property
    def min_tile_num(self):
        min_tile = self._args.get('--min_tile')
        # extra logic in case min tile needs to be zero
        if min_tile is False:
            min_tile = 1
        return 2100 + int(min_tile)

    @property
    def min_hits(self):
        return int(self._args.get('--min_hits') or 15)

    @property
    def objective(self):
        return int(self._args.get('--objective') or 60)

    @property
    def chip_id(self):
        return self._args.get('--chip_id') or self._presumptive_chip_id

    @property
    def rotation_estimate(self):
        return float(self._args.get('--rotation_estimate') or 180.0)

    @property
    def snr_threshold(self):
        return float(self._args.get('--snr_threshold') or 1.2)

    def _make_filepath(self, filename):
        return '{base_directory}{sep}{chip_id}{sep}{filename}'.format(base_directory=self._base_directory,
                                                                      chip_id=self.chip_id,
                                                                      sep=os.path.sep,
                                                                      filename=filename)
