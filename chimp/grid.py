import re
import numpy as np


class Image(np.ndarray):
    """
    Holds the raw pixel data of an image and provides access to some metadata.

    """
    def __new__(cls, array, row, column):
        return np.asarray(array).view(cls)

    def __init__(self, array, row, column):
        self.row = row
        self.column = column

    @property
    def index(self):
        return "%.3d_%.3d" % (self.row, self.column)

    def __array_wrap__(self, obj, *_):
        if len(obj.shape) == 0:
            return obj[()]
        else:
            return obj


class GridImages(object):
    def __init__(self, h5, channel):
        """
        Since some ND2s were created where multiple channels had the same name, we can't always use the channel name,
        though we will be able to going forward now that we saw that that was happening.

        In the near future, we will also have data sets where we take multiple images in the same location and in the
        same channel, with different exposures. We will then need to add them together, dynamically choosing which ones
        to add in order to avoid saturation.

        """
        self._h5 = h5
        self._rows = None
        self._columns = None
        self._channel = channel
        self._parse_grid()

    def __iter__(self):
        for image in self.bounded_iter(0, self._width):
            yield image

    def _parse_grid(self):
        regex = re.compile('''^\(Major, minor\) = \((?P<column>\d+), (?P<row>\d+)\)$''')
        max_row = 0
        max_column = 0
        for key in self._h5[self._channel].keys():
            match = regex.search(key)
            if match:
                column, row = match.group('column'), match.group('row')
                max_column = max(max_column, int(column))
                max_row = max(max_row, int(row))
        self._height = max_row + 1
        self._width = max_column + 1

    def bounded_iter(self, min_column, max_column):
        """
        Iterates over all images between two columns (inclusive)

        """
        for column in range(min_column, max_column):
            for row in range(self._height):
                image = self.get(row, column)
                if image is not None:
                    yield image

    def left_iter(self):
        return self.bounded_iter(0, self._width)

    def right_iter(self):
        for column in reversed(range(self._width)):
            for row in reversed(range(self._height)):
                image = self.get(row, column)
                if image is not None:
                    yield image

    def get(self, row, column):
        try:
            raw_array = self._h5[self._channel]['(Major, minor) = (%d, %d)' % (column, row)].value
            return Image(raw_array, row, column)
        except (KeyError, IndexError, AttributeError):
            return None
