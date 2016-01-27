import os


class FileStructure(object):
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
