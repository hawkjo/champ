import numpy as np


class XYZFile(object):
    def __init__(self, image):
        self._image = image

    def __str__(self):
        column_width = len(str(self._image.shape[0]))
        row_width = len(str(self._image.shape[1]))
        intensity_width = len(str(np.max(self._image)))
        line = "{column: >%s} {row: >%s} {intensity: >%s}" % (column_width, row_width, intensity_width)
        return "\n".join(line.format(row=row,
                                     column=column,
                                     intensity=intensity) for row, column, intensity in self._as_vector())

    def _as_vector(self):
        for (row, column), intensity in np.ndenumerate(self._image):
            yield row, column, intensity
