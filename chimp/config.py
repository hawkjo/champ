import os


class Experiment(object):
    def __init__(self, project_name):
        self.project_name = project_name

    @property
    def figure_directory(self):
        return 'figs'

    @property
    def data_directory(self):
        return 'data'

    @property
    def fourier_data_directory(self):
        return os.path.join(self.data_directory, 'from_fourierseq')

    @property
    def results_directory(self):
        return 'results'

    def get_sexcat_path(self, nd2_name, image_index):
        return os.path.join(nd2_name, '%d.cat' % image_index)


class AlignmentParameters(object):
    """ Parses user-provided alignment parameters and provides a default in case no value was given. """
    def __init__(self, command_line_args):
        self._args = command_line_args

    @property
    def aligning_read_names_filepath(self):
        return 'phiX_mappings/phiX_read_names.txt'

    @property
    def all_read_names_filepath(self):
        return 'read_names/all_read_names.txt'

    @property
    def fq_w_est(self):
        return 935.0

    @property
    def min_hits(self):
        return 15

    @property
    def rotation_estimate(self):
        return 180.0

    @property
    def snr_threshold(self):
        return 1.2
