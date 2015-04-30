import numpy as np
import matplotlib.pyplot as plt
from image_processing_tools import max_2d_idx


def unique_rows(a):
    """Returns unique rows of a matrix."""
    assert a.ndim == 2, 'unique_rows method works only for 2d arrays.'
    b = a.view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    return np.unique(b).view(a.dtype).reshape(-1, a.shape[1])


class ImagePointsCorrelation(object):
    """A class for computing the cross-correlation between an image and a field of points.
    Useful when the number of points is small relative to the number of pixels required to
    represent them as an image."""
    def __init__(self, im, points):
        self.im = im
        self.points = points
        self.find_points_bounding_box()
        self.find_point_pixels()
        self.correlate()

    def correlate(self):
        """Computes the raw cross-correlation between the image and the point-field."""
        # The raw cross-correlation is a bit of a strange object. Consider an image and a
        # point-field which fit in the following bounding boxes:
        #
        #   Image                           Bounding box of point field
        #   --------------------            ----------
        #   |                  |            |        |
        #   |                  |            |        |
        #   |                  |            ----------
        #   |                  |
        #   |                  |
        #   --------------------
        #
        # Then the raw cross-correlation will look like this:
        #
        #   ------------------------------
        #   |         *                  |
        #   |         *                  |
        #   |         *                  |
        #   |*********O******************|
        #   |         *                  |
        #   |         *                  |
        #   |         *                  |
        #   |         *                  |
        #   |         *                  |
        #   ------------------------------
        #
        # where the 'O' represents the origin of the cross-correlation, even though it is not in
        # element [0, 0] of the matrix. Specifically, it is in the row and column equal to the
        # dimensions of the bounding box of the point field.
        #
        # Through proper mappings, this gives a natural way to represent correlation offsets both
        # positive and negative. Note that those mappings are 
        #
        #   1) Correct for offset of the origin
        #   2) Correct for offset to lower-left corner of the bounding box of the point field.

        imshape = self.im.shape
        n_rows = imshape[0] + self.point_pixels_im_shape[0]
        n_cols = imshape[1] + self.point_pixels_im_shape[1]

        self.raw_correlation = np.zeros((n_rows, n_cols))
        self.origin = self.point_pixels_im_shape
        for point_px in self.point_pixels:
            cxy = self.origin - point_px
            self.raw_correlation[cxy[0]:cxy[0]+imshape[0], cxy[1]:cxy[1]+imshape[1]] += self.im
        self.pixel_tr = max_2d_idx(self.raw_correlation) - self.origin
        self.point_tr = self.pixel_tr - self.q1_offset

    def find_points_bounding_box(self, force=False):
        if hasattr(self, 'y_width') and not force:
            return
        self.x_min, self.y_min = self.points.min(axis=0)
        self.x_max, self.y_max = self.points.max(axis=0)
        self.x_width = self.x_max - self.x_min
        self.y_width = self.y_max - self.y_min
        self.q1_offset = np.array([-self.x_min, -self.y_min])

    def find_point_pixels(self):
        """Moves points to first quadrant/origin, converts to integers, and uniques list."""
        self.point_pixels = unique_rows(
                (self.points + np.tile(self.q1_offset, (self.points.shape[0], 1))).astype(np.int)
                )
        self.point_pixels_im_shape = self.point_pixels.max(axis=0) + [1, 1]

    def find_translated_points(self):
        self.translated_points = self.points + np.tile(self.point_tr, (self.points.shape[0], 1))

    def overlay_points_on_image(self, ax=None):
        assert hasattr(self, 'raw_correlation') and self.raw_correlation.size, \
                'Correlation not computed.'
        if not hasattr(self, 'translated_points'):
            self.find_translated_points()

        if ax is None:
            fig, ax = plt.subplots()
        ax.imshow(self.im, cmap=plt.get_cmap('Blues'))
        ax.plot(self.translated_points[:, 1], self.translated_points[:, 0],
                'r.',
                alpha=0.4,
                linewidth=0)

    def imshow_correlation(self, ax=None):
        assert hasattr(self, 'raw_correlation') and self.raw_correlation.size, \
                'Correlation not computed.'
        if not hasattr(self, 'translated_points'):
            self.find_translated_points()

        if ax is None:
            fig, ax = plt.subplots()
        ax.imshow(self.raw_correlation)
        ax.plot([self.origin[0]], [self.origin[1]], 'w*')
