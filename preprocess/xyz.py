import numpy as np


class XYZFile(object):
    def __init__(self, image):
        self._image = image

    def __str__(self):
        return "\n".join("{column} {row} {intensity}".format(row=row, column=column, intensity=intensity) for row, column, intensity in self._as_vector())

    def _as_vector(self):
        for (row, column), intensity in np.ndenumerate(self._image):
            yield row, column, intensity
