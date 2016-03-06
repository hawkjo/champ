import numpy as np


class SextractorPoint(object):
    __slots__ = ['column', 'row', 'flux', 'flux_err', 'flags', 'width', 'height', 'theta']

    def __init__(self, column, row, flux, flux_err, flags, width, height, theta):
        # Sextractor coordinates are 1-based
        self.column = column - 1
        self.row = row - 1
        self.flux = flux
        self.flux_err = flux_err
        self.flags = flags
        self.width = width
        self.height = height
        self.theta = theta


class Sextraction(object):
    def __init__(self, lines):
        self.points = []
        for line in filter(lambda x: not x.startswith('#'), lines):
            self.points.append(SextractorPoint(*line.strip().split()))

    @property
    def rcs(self):
        return np.array([(pt.row, pt.column) for pt in self.points])
