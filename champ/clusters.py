import numpy as np


def parse_sextractor(line):
    column, row, _, _, _, _, _, _ = map(float, line.strip().split())
    # Sextractor coordinates are 1-based
    return row - 1.0, column - 1.0


def parse_cluster(line):
    row, column = map(float, line.strip().split())
    return row, column


class ClusterPoint(object):
    __slots__ = ('row', 'column')

    def __init__(self, row, column):
        self.row = float(row)
        self.column = float(column)


class Clusters(object):
    def __init__(self, lines, parser_name):
        parsers = {'sextractor': parse_sextractor,
                   'otsu': parse_cluster}
        parse = parsers[parser_name]
        self.points = []
        for line in lines:
            row, column = parse(line)
            self.points.append(ClusterPoint(row, column))
        self.point_rcs = np.array([(pt.row, pt.column) for pt in self.points])

    def rs(self):
        return np.array([pt.row for pt in self.points])

    def cs(self):
        return np.array([pt.column for pt in self.points])
