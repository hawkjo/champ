import os
import logging


class CommandLineArguments(object):
    """
    Wraps the raw arguments provided by docopt.

    """
    def __init__(self, arguments, current_directory):
        self._arguments = arguments
        self._current_directory = current_directory

    @property
    def log_level(self):
        log_level = {0: logging.FATAL,
                     1: logging.WARN,
                     2: logging.INFO,
                     3: logging.DEBUG}
        # default to silent if the user supplies no verbosity setting
        return log_level.get(self._arguments['-v'], logging.FATAL)

    @property
    def alignment_channel(self):
        return self._arguments['ALIGNMENT_CHANNEL']

    @property
    def image_directory(self):
        return self._arguments['IMAGE_DIRECTORY']

    @property
    def project_name(self):
        return self._arguments['PROJECT_NAME']

    @property
    def fastq_directory(self):
        return self._arguments['FASTQ_DIRECTORY']

    @property
    def hdf5_file_path(self):
        return self._arguments['HDF5_FILE_PATH']

    @property
    def tif_directories(self):
        return self._arguments['TIF_DIRECTORIES']

    @property
    def bamfiles(self):
        return self._arguments['PATHS_TO_BAMFILES']

    @property
    def command(self):
        for possible_command in ('align',
                                 'preprocess',
                                 'map',
                                 'convert'):
            if self._arguments[possible_command]:
                return possible_command

    @property
    def snr_threshold(self):
        return float(self._arguments['--snr-threshold'])

    @property
    def min_hits(self):
        return int(self._arguments['--min-hits'])

    @property
    def flipud(self):
        # flip images across the horizontal axis
        return self._arguments['--flipud']

    @property
    def fliplr(self):
        # flip images across the horizontal axis
        return self._arguments['--fliplr']


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
    def alignment_file(self):
        return os.path.join(self.results_directory, 'align')

    @property
    def results_directory(self):
        return 'results'


class AlignmentParameters(object):
    """ Parses user-provided alignment parameters and provides a default in case no value was given. """
    def __init__(self, command_line_args):
        self._args = command_line_args

    @property
    def aligning_read_names_filepath(self):
        return 'mapped_reads/phix'

    @property
    def all_read_names_filepath(self):
        return 'mapped_reads/unclassified'

    @property
    def fastq_tile_width_estimate(self):
        # width of a tile of Illumina data, in microns
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
