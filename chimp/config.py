import os
import logging
from chip import Miseq, Hiseq


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
    def second_channel(self):
        return self._arguments['SECOND_CHANNEL_NAME']

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
    def microns_per_pixel(self):
        return float(self._arguments['MICRONS_PER_PIXEL'])

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
    def chip(self):
        chip = self._arguments.get('CHIP_TYPE', 'miseq')
        chips = {'miseq': Miseq,
                 'hiseq': Hiseq}
        return chips[chip](self._arguments['--ports-on-left'])

    @property
    def snr_threshold(self):
        return float(self._arguments.get('SNR', 1.2))

    @property
    def min_hits(self):
        return int(self._arguments.get('MIN_HITS', 15))

    @property
    def make_pdfs(self):
        return self._arguments['--make-pdfs']

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
        #
        return 935.0

    @property
    def min_hits(self):
        return self._args.min_hits

    @property
    def rotation_estimate(self):
        return 180.0

    @property
    def snr_threshold(self):
        return self._args.snr_threshold
