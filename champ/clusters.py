import numpy as np


class ClusterPoint(object):
    __slots__ = ('r', 'c')

    def __init__(self, line):
        self.r, self.c = map(float, line.strip().split())


class Clusters(object):
    def __init__(self, lines):
        self.points = []
        for line in lines:
            self.points.append(ClusterPoint(line))
        self.point_rcs = np.array([(pt.r, pt.c) for pt in self.points])

    def rs(self):
        return np.array([pt.r for pt in self.points])

    def cs(self):
        return np.array([pt.c for pt in self.points])
