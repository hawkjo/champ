from copy import deepcopy
import numpy as np
import misc
import logging
from scipy import ndimage

log = logging.getLogger(__name__)


class FastqTileRCs(object):
    """A class for fastq tile coordinates."""
    def __init__(self, key, read_names, microns_per_pixel):
        self.key = key
        self.microns_per_pixel = microns_per_pixel
        self.read_names = read_names
        self.rcs = np.array([map(int, name.split(':')[-2:]) for name in self.read_names])

    def set_fastq_image_data(self, offset, scale, scaled_dims, width):
        self.offset = offset
        self.scale = scale
        self.image_shape = scaled_dims
        self.width = width  # width in um
        self.mapped_rcs = scale * (self.rcs + np.tile(offset, (self.rcs.shape[0], 1)))
        self.rotation_degrees = 0

    def rotate_data(self, degrees):
        self.rotation_degrees += degrees
        self.mapped_rcs = np.dot(self.mapped_rcs, misc.right_rotation_matrix(degrees, degrees=True))
        self.mapped_rcs -= np.tile(self.mapped_rcs.min(axis=0), (self.mapped_rcs.shape[0], 1))
        self.image_shape = self.mapped_rcs.max(axis=0) + 1
        return self.image_shape

    def image(self):
        image = np.zeros(self.image_shape.astype(np.int))
        image[self.mapped_rcs.astype(np.int)[:, 0], self.mapped_rcs.astype(np.int)[:, 1]] = 1
        sigma = 0.25 / self.microns_per_pixel  # Clusters have stdev ~= 0.25 um
        image = ndimage.gaussian_filter(image, sigma)
        return image

    def fft_align_with_im(self, image_data):
        # Make the ffts
        fq_image = self.image()
        padded_fq_im = misc.pad_to_size(fq_image, image_data.fft.shape)
        fq_im_fft = np.fft.fft2(padded_fq_im)
        # Align
        im_data_fft = image_data.fft
        cross_corr = abs(np.fft.ifft2(np.conj(fq_im_fft) * im_data_fft))
        max_corr = cross_corr.max()
        max_idx = misc.max_2d_idx(cross_corr)
        align_tr = np.array(max_idx) - fq_image.shape
        return max_corr, align_tr

    def set_aligned_rcs(self, align_tr):
        log.debug("ALIGN TR: %s" % align_tr)
        """Returns aligned rcs. Only works when image need not be flipped or rotated."""
        self.aligned_rcs = deepcopy(self.mapped_rcs)
        self.aligned_rcs -= np.tile(self.aligned_rcs.min(axis=0), (self.aligned_rcs.shape[0], 1))
        self.aligned_rcs += np.tile(align_tr, (self.aligned_rcs.shape[0], 1))

    def set_aligned_rcs_given_transform(self, lbda, theta, offset):
        """Performs transform calculated in FastqImageCorrelator.least_squares_mapping."""
        A = np.zeros((2*len(self.rcs), 4))
        for i, pt in enumerate(self.rcs):
            xir, yir = pt
            A[2*i, :] = [xir, -yir, 1, 0]
            A[2*i+1, :] = [yir,  xir, 0, 1]

        x = np.array([lbda * np.cos(theta),
                      lbda * np.sin(theta),
                      offset[0],
                      offset[1]])

        # First update w since it depends on previous scale setting
        self.width = lbda * float(self.width) / self.scale
        self.scale = lbda
        self.rotation = theta
        self.rotation_degrees = theta * 180.0 / np.pi
        self.offset = offset
        self.aligned_rcs = np.dot(A, x).reshape((len(self.rcs), 2))

    def set_correlation(self, im):
        """Sets alignment correlation. Only works when image need not be flipped or rotated."""
        self.best_max_corr = sum(im[int(x), int(y)] for x, y in self.aligned_rcs
                                 if 0.0 <= x < im.shape[0] and 0.0 <= y < im.shape[1])

    def set_snr_with_control_corr(self, control_corr):
        self.snr = self.best_max_corr / control_corr
