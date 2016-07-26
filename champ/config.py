import os
import logging
from chip import Miseq, Hiseq, load


class CommandLineArguments(object):
    """
    Wraps the raw arguments provided by docopt.

    """
    def __init__(self, arguments, current_directory):
        self._arguments = arguments
        self._current_directory = current_directory

    @property
    def log_level(self):
        log_level = {0: logging.ERROR,
                     1: logging.WARN,
                     2: logging.INFO,
                     3: logging.DEBUG}
        # default to silent if the user supplies no verbosity setting
        return log_level.get(self._arguments['-v'], logging.ERROR)

    @property
    def alignment_channel(self):
        return self._arguments['ALIGNMENT_CHANNEL']

    @property
    def nonneg_lda_weights_path(self):
        return 'bLDA_coef_nonneg.txt'

    @property
    def image_directory(self):
        return self._arguments['IMAGE_DIRECTORY']

    @property
    def fastq_directory(self):
        return self._arguments['FASTQ_DIRECTORY']

    @property
    def mapped_reads(self):
        return self._arguments['MAPPED_READS']

    @property
    def hdf5_file_path(self):
        return self._arguments['HDF5_FILE_PATH']

    @property
    def microns_per_pixel(self):
        mpp = self._arguments.get('--microns-per-pixel') or 0.2666666666666666666
        return float(mpp)

    @property
    def output_directory(self):
        return self._arguments['OUTPUT_DIRECTORY']

    @property
    def bamfiles(self):
        # TODO: Auto-expand user directories or convert these to absolute paths or something
        return self._arguments['PATHS_TO_BAMFILES']

    @property
    def command(self):
        for possible_command in ('map',
                                 'init',
                                 'align',
                                 'kd',
                                 'info'):
            if self._arguments.get(possible_command):
                return possible_command

    @property
    def chip(self):
        chip = load(self._arguments.get('CHIP_TYPE', 'miseq'))
        return chip(self._arguments['--ports-on-right'])

    @property
    def phix_only(self):
        return self._arguments['--phix-only']

    @property
    def ports_on_right(self):
        return self._arguments['--ports-on-right']

    @property
    def chip_name(self):
        return self._arguments['CHIP_NAME']

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

    @property
    def force(self):
        return self._arguments['--force']


class Experiment(object):
    def __init__(self, image_data_directory):
        self._image_data_directory = image_data_directory

    @property
    def figure_directory(self):
        return os.path.join(self._image_data_directory, 'figs')

    @property
    def data_directory(self):
        return self._image_data_directory

    @property
    def results_directory(self):
        return os.path.join(self._image_data_directory, 'results')

    @property
    def intensity_directory(self):
        return os.path.join(self.figure_directory, 'intensity')


class AlignmentParameters(object):
    """ Parses user-provided alignment parameters and provides a default in case no value was given. """
    def __init__(self, clargs, mapped_reads):
        self._clargs = clargs
        self._mapped_reads = mapped_reads

    @property
    def aligning_read_names_filepath(self):
        return os.path.join(self._mapped_reads, 'phix')

    @property
    def all_read_names_filepath(self):
        return os.path.join(self._mapped_reads, 'unclassified')

    @property
    def fastq_tile_width_estimate(self):
        # width of a tile of Illumina data, in microns
        return 935.0

    @property
    def min_hits(self):
        return self._clargs.min_hits

    @property
    def rotation_estimate(self):
        return 180.0

    @property
    def snr_threshold(self):
        return self._clargs.snr
