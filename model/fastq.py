class FastqRead(object):
    __slots__ = ['name', 'sequence', '_column', '_row', '_lane', '_tile']

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self._lane, self._tile, self._column, self._row = name.rsplit(':')[-4:]

    @property
    def tile(self):
        return tuple(map(int, (self._lane, self._tile[3:5])))

    @property
    def row(self):
        return int(self._row)

    @property
    def column(self):
        return int(self._column)
