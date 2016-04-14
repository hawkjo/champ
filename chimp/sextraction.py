import numpy as np


class SextractorPoint(object):
    def __init__(self, line):
        self.c, self.r, self.flux, self.flux_err, self.flags, \
                self.width, self.height, self.theta \
                = map(float, line.strip().split())
        # Sextractor coordinates are 1-based
        self.r -= 1
        self.c -= 1


class Sextraction(object):
    def __init__(self, lines):
        self.points = []
        for line in lines:
            if line.startswith('#'):
                continue
            self.points.append(SextractorPoint(line))
        self.point_rcs = np.array([(pt.r, pt.c) for pt in self.points])

    def rs(self):
        return np.array([pt.r for pt in self.points])

    def cs(self):
        return np.array([pt.c for pt in self.points])
