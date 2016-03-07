import logging
import os


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
    def out_directory(self):
        return self._arguments['-o'] if self._arguments['-o'] else os.path.join(self._current_directory,
                                                                                'classified_reads')

    @property
    def fastq_directory(self):
        return self._arguments['FASTQ_DIRECTORY']

    @property
    def bamfiles(self):
        return self._arguments['PATHS_TO_BAMFILES']

    @property
    def command(self):
        for possible_command in ('align',
                                 'bowtie',
                                 'readsort'):
            if self._arguments[possible_command]:
                return possible_command

    @property
    def min_hits(self):
        return int(self._arguments['--min_hits'])

    @property
    def tile_width_estimate(self):
        return int(self._arguments['--tile_width_estimate'])

    @property
    def rotation_estimate(self):
        return float(self._arguments['--rotation_estimate'])

    @property
    def index_offset(self):
        # TODO: This might be unnecessary
        return int(self._arguments['--index_offset'])
