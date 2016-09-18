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
    def cluster_strategy(self):
        return self._arguments.get('CLUSTER_STRATEGY', 'sextractor')

    @property
    def nonneg_lda_weights_path(self):
        return self._arguments['LDA_WEIGHTS']

    @property
    def perfect_target_name(self):
        return self._arguments.get('PERFECT_TARGET_NAME', False)

    @property
    def alternate_fiducial_reads(self):
        # sometimes you don't want to use phix for end tile finding and initial alignment
        return self._arguments.get('ALTERNATE_FIDUCIAL_READS', False)

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
        return self._arguments['READ_NAMES_DIRECTORY']

    @property
    def parsed_reads(self):
        # TODO: This needs to be merged into mapped_reads
        return self._arguments['READ_NAMES_DIRECTORY']

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
        # We have to do this weird loop to deal with the way docopt stores the command name
        for possible_command in ('map',
                                 'init',
                                 'align',
                                 'info'):
            if self._arguments.get(possible_command):
                return possible_command

    @property
    def chip(self):
        chip = load(self._arguments.get('--chip') or 'miseq')
        return chip(self._arguments['--ports-on-right'])

    @property
    def fiducial_only(self):
        return self._arguments['--fiducial-only']

    @property
    def ports_on_right(self):
        return self._arguments['--ports-on-right']

    @property
    def chip_name(self):
        return self._arguments['CHIP_NAME']

    @property
    def min_hits(self):
        return int(self._arguments.get('MIN_HITS') or 15)

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
    def alternate_good_target_reads_filename(self):
        return self._arguments.get('ALTERNATE_GOOD_READS')

    @property
    def alternate_perfect_target_reads_filename(self):
        return self._arguments.get('ALTERNATE_PERFECT_READS')

    @property
    def force(self):
        return self._arguments['--force']

    @property
    def snr(self):
        return float(self._arguments.get('SNR') or 1.2)


class PathInfo(object):
    """ Parses user-provided alignment parameters and provides a default in case no value was given. """
    def __init__(self, image_directory, mapped_reads, perfect_target_name, alternate_fiducial_reads=None,
                 alternate_perfect_reads_filename=None, alternate_good_reads_filename=None):
        self._image_directory = image_directory
        self._mapped_reads = mapped_reads
        self._alternate_fiducial_reads = alternate_fiducial_reads
        self._perfect_target_name = perfect_target_name
        self._alternate_perfect_reads_filename = alternate_perfect_reads_filename
        self._alternate_good_reads_filename = alternate_good_reads_filename

    @property
    def figure_directory(self):
        return os.path.join(self._image_directory, 'figs')

    @property
    def results_directory(self):
        return os.path.join(self._image_directory, 'results')

    @property
    def aligning_read_names_filepath(self):
        if self._alternate_fiducial_reads:
            return os.path.join(self._mapped_reads, self._alternate_fiducial_reads)
        return os.path.join(self._mapped_reads, 'phix')

    @property
    def all_read_names_filepath(self):
        return os.path.join(self._mapped_reads, 'all_read_names.txt')

    @property
    def on_target_read_names(self):
        if self._alternate_good_reads_filename:
            return os.path.join(self._mapped_reads, self._alternate_good_reads_filename)
        if not self._perfect_target_name:
            raise ValueError("This experiment did not have a perfect target set!")
        return os.path.join(self._mapped_reads, 'target_{}_read_names.txt'.format(self._perfect_target_name.lower()))

    @property
    def perfect_read_names(self):
        if self._alternate_perfect_reads_filename:
            return os.path.join(self._mapped_reads, self._alternate_perfect_reads_filename)
        if not self._perfect_target_name:
            raise ValueError("This experiment did not have a perfect target set!")
        return os.path.join(self._mapped_reads, 'perfect_target_{}_read_names.txt'.format(self._perfect_target_name.lower()))
