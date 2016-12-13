import numpy as np


class ClusterPoint(object):
    __slots__ = ('r', 'c')

    def __init__(self, line):
        self.r, self.c = map(float, line.strip().split())


class SextractorPoint(object):
    def __init__(self, line):
        self.c, self.r, self.flux, self.flux_err, self.flags, \
                self.width, self.height, self.theta \
                = map(float, line.strip().split())
        # Sextractor coordinates are 1-based
        self.r -= 1
        self.c -= 1


class Clusters(object):
    def __init__(self, lines, cluster_strategy):
        # point objects are essentially parsers that wrap the underlying coordinates of clusters
        point_objects = {'otsu': ClusterPoint,
                         'se': SextractorPoint}
        Point = point_objects[cluster_strategy]
        self.points = []
        for line in lines:
            if line.startswith("#"):
                continue
            self.points.append(Point(line))
        self.point_rcs = np.array([(pt.r, pt.c) for pt in self.points])

    def rs(self):
        return np.array([pt.r for pt in self.points])

    def cs(self):
        return np.array([pt.c for pt in self.points])
