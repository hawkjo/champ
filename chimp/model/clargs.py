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
    def project_name(self):
        return self._arguments['PROJECT_NAME']

    @property
    def fastq_directory(self):
        return self._arguments['FASTQ_DIRECTORY']

    @property
    def bamfiles(self):
        return self._arguments['PATHS_TO_BAMFILES']

    @property
    def alignment_channel(self):
        return self._arguments['--alignment-channel']

    @property
    def command(self):
        for possible_command in ('align',
                                 'preprocess'):
            if self._arguments[possible_command]:
                return possible_command

    @property
    def snr_threshold(self):
        return float(self._arguments['--snr-threshold'])

    @property
    def min_hits(self):
        return int(self._arguments['--min-hits'])
