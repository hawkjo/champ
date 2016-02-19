class FastqRead(object):
    __slots__ = ['_record', '_lane', '_tile', '_column', '_row']

    def __init__(self, record):
        self._record = record
        self._lane, self._tile, self._column, self._row = map(int, record.name.rsplit(':')[-4:])

    def __repr__(self):
        return self._record.format("fastq")

    @property
    def name(self):
        return self._record.name

    @property
    def tile(self):
        return self._lane, self._tile

    @property
    def row(self):
        return self._row

    @property
    def column(self):
        return self._column
