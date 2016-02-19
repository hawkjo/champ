class FastqFiles(object):
    """ Sorts compressed FastQ files provided to us from the Illumina sequencer. """
    def __init__(self, filenames):
        self._filenames = [f for f in self._filter_names(filenames)]

    @property
    def paired_files(self):
        for f1, f2 in self._sort_filenames(paired=True):
            yield f1, f2

    @property
    def single_files(self):
        for f in self._sort_filenames(paired=False):
            yield f

    def _filter_names(self, filenames):
        for filename in filenames:
            if not filename.endswith('fastq.gz'):
                continue
            for text in ('_I1_', '_I2_', '_I1.', '_I2.'):
                if text in filename:
                    continue
            yield filename

    def _sort_filenames(self, paired=True):
        for filename in self._filenames:
            if '_R1_' in filename or '_R1.' in filename:
                pair = filename.replace('_R1_', '_R2_').replace('_R1.', '_R2.')
                if paired and pair in self._filenames:
                    yield filename, pair
                elif not paired and pair not in self._filenames:
                    yield filename
