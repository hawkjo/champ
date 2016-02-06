import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


class SextractorPoint(object):
    __slots__ = ['c', 'r', 'flux', 'flux_err', 'flags', 'width', 'height', 'theta']

    def __init__(self, line):
        self.c, self.r, self.flux, self.flux_err, self.flags, \
                self.width, self.height, self.theta \
                = map(float, line.strip().split())
        # Sextractor coordinates are 1-based
        self.r -= 1
        self.c -= 1


class Sextraction(object):
    def __init__(self, fpath):
        self.points = []
        for line in open(fpath):
            if line.startswith('#'):
                continue
            self.points.append(SextractorPoint(line))
        self.point_rcs = np.array([(pt.r, pt.c) for pt in self.points])

    def rs(self):
        return np.array([pt.r for pt in self.points])

    def cs(self):
        return np.array([pt.c for pt in self.points])

    def plot_ellipses(self, ax=None, alpha=1.0, color=(1, 0, 0)):
        if ax is None:
            fig, ax = plt.subplots()
        ells = [Ellipse(xy=(pt.c, pt.r), width=pt.width, height=pt.height, angle=pt.theta) for pt in self.points]
        for e in ells:
            ax.add_artist(e)
            e.set_alpha(alpha)
            e.set_facecolor(color)
