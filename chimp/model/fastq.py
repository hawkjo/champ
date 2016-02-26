class FastqRead(object):
    """ Wraps the raw data about a single DNA read that we receive from Illumina. """
    __slots__ = ('label', '_name', '_seq', '_lane', '_tile', '_column', '_row')

    def __init__(self, record):
        self._name = record.name
        self._seq = record.seq
        self._lane, self._tile, self._column, self._row = map(int, record.name.rsplit(':')[-4:])
        self.label = None

    @property
    def name(self):
        return self._name

    @property
    def sequence(self):
        return self._seq

    @property
    def region(self):
        return self._lane, self._tile

    @property
    def row(self):
        return self._row

    @property
    def column(self):
        return self._column


class FastqFiles(object):
    """ Sorts compressed FastQ files provided to us from the Illumina sequencer. """
    def __init__(self, filenames):
        self._filenames = filenames
        
    def __iter__(self):
        for f in self._filenames:
            yield f

    @property
    def paired_files(self):
        for f1, f2 in self._sort_filenames(paired=True):
            yield f1, f2

    @property
    def single_files(self):
        for f in self._sort_filenames(paired=False):
            yield f

    def _filter_names(self, data):
        for filename in data:
            if not filename.endswith('fastq.gz'):
                continue
            if '_I1_' in filename or '_I2_' in filename or '_I1.' in filename or '_I2.' in filename:
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
