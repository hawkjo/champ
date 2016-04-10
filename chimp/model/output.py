from abc import abstractmethod
from chimp.model.constants import FASTQ_TILE_WIDTH
import math


class ResultFile(object):
    def __init__(self, index, filename):
        self._index = int(index)
        self._filename = filename
        self._lines = set()

    @property
    def filename(self):
        return "%d_%s.txt" % (self._index, self._filename)

    @abstractmethod
    def __str__(self):
        raise NotImplementedError


class Intensities(ResultFile):
    def __init__(self, index):
        super(Intensities, self).__init__(index, 'intensities')


class Stats(ResultFile):
    def __init__(self, index, project_name, objective, lane, side, results):
        super(Stats, self).__init__(index, 'stats')
        self._project_name = project_name
        self._objective = objective
        self._lane = lane
        self._side = side
        self._results = results

    def _tile_keys(self):
        """ Converts tile information into legacy tile names. """
        for result in self._results:
            yield 'lane%dtile%d1%.2d' % (result.lane, result.side, result.tile_number)

    def __str__(self):
        return '\n'.join(('Image:                 %s' % self._index,
                          'Objective:             %d' % self._objective,
                          'Project Name:          %s' % self._project_name,
                          'Tile:                  %s' % ','.join(self._tile_keys()),
                          'Rotation (deg):        %s' % ','.join('%.4f' % (180.0 + result.theta * 180.0 / math.pi) for result in self._results),
                          'Tile width (um):       %s' % ','.join('%.4f' % FASTQ_TILE_WIDTH for _ in self._results),
                          'Scaling (px/fqu):      %s' % ','.join('%.7f' % result.scale for result in self._results),
                          'RC Offset (px):        %s' % ','.join('(%.4f,%.4f)' % tuple(result.offset) for result in self._results),
                          'Correlation:           %s' % ','.join('%.2f' % result.correlation for result in self._results),
                          'Non-mutual hits:       %d' % (sum([len(result.hits.non_mutual) for result in self._results])),
                          'Bad mutual hits:       %d' % (sum([len(result.hits.bad_mutual) for result in self._results])),
                          'Good mutual hits:      %d' % (sum([len(result.hits.good_mutual) for result in self._results])),
                          'Exclusive hits:        %d' % (sum([len(result.hits.exclusive) for result in self._results])),
                          'Sextractor Ellipses:   %d' % (sum([len(result.sexcat_rcs) for result in self._results])),
                          'Fastq Points:          %d' % sum(result.fastq_point_count for result in self._results)))


class AllReadRCs(ResultFile):
    """
    Holds and formats results that go into the *_all_reads_rcs.txt file

    """
    def __init__(self, index):
        super(AllReadRCs, self).__init__(index, 'all_read_rcs')

    def add_rcs(self, original_rcs, aligned_rcs, fastq_read_rcs):
        for (original_r, original_c), (aligned_r, aligned_c) in zip(original_rcs, aligned_rcs):
            name = fastq_read_rcs.get((original_r, original_c))
            if name is not None:
                self._lines.add((name, aligned_r, aligned_c))

    def __str__(self):
        return '\n'.join(["%s\t%.6f\t%.6f" % (name, r, c) for name, r, c in self._lines])
