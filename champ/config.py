import logging
import os

from chip import load


class CommandLineArguments(object):
    """
    Wraps the raw arguments provided by docopt.

    """
    def __init__(self, arguments, current_directory):
        self._arguments = arguments
        self._current_directory = current_directory

    @property
    def target_data_file(self):
        return self._arguments['TARGET_DATA_FILE']

    @property
    def nonneg_lda_weights_path(self):
        return self._arguments['LDA_WEIGHTS']

    @property
    def target_label(self):
        return self._arguments['TARGET_LABEL']

    @property
    def off_target_label(self):
        return self._arguments['OFF_TARGET_LABEL']

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
    def image_directory(self):
        return self._arguments['IMAGE_DIRECTORY']

    @property
    def fastq_directory(self):
        return self._arguments['FASTQ_DIRECTORY']

    @property
    def mapped_reads(self):
        return self._arguments['MAPPED_READS']

    @property
    def parsed_reads(self):
        return self._arguments['PARSED_READS']

    @property
    def microns_per_pixel(self):
        return float(self._arguments.get('--microns-per-pixel') or 0.2666666666666666666)

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
                                 'preprocess',
                                 'info'):
            if self._arguments.get(possible_command):
                return possible_command

    @property
    def chip(self):
        chip = load(self._arguments.get('--chip') or 'miseq')
        return chip(self._arguments['--ports-on-right'])

    @property
    def ports_on_right(self):
        return self._arguments['--ports-on-right']

    @property
    def chip_name(self):
        return self._arguments['CHIP_NAME']

    @property
    def min_hits(self):
        return int(self._arguments.get('--min-hits') or 15)

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

    @property
    def snr(self):
        return float(self._arguments.get('-snr') or 1.2)


class OutputParameters(object):
    # TODO: The name of this is bad, it's input and output
    """ Parses user-provided alignment parameters and provides a default in case no value was given. """
    def __init__(self, image_directory, mapped_reads):
        self._image_directory = image_directory
        self._mapped_reads = mapped_reads

    @property
    def figure_directory(self):
        return os.path.join(self._image_directory, 'figs')

    @property
    def results_directory(self):
        return os.path.join(self._image_directory, 'results')

    @property
    def aligning_read_names_filepath(self):
        return os.path.join(self._mapped_reads, 'phix')

    @property
    def all_read_names_filepath(self):
        return os.path.join(self._mapped_reads, 'unclassified')
