from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.optimize import minimize
from scipy.ndimage.filters import gaussian_filter
import misc 


class FastqTileRCs(object):
    """A class for fastq tile coordinates."""
    def __init__(self, key, read_names):
        self.key = key
        self.read_names = read_names
        self.rcs = np.array([map(int, name.split(':')[-2:]) for name in self.read_names])

    def set_fastq_image_data(self, offset, scale, scaled_dims, w, um_per_pixel, force=False, verbose=True):
        self.offset = offset
        self.scale = scale
        self.image_shape = map(int, scaled_dims)
        self.w = w  # width in um
        self.um_per_pixel = um_per_pixel
        self.mapped_rcs = scale * (self.rcs + np.tile(offset, (self.rcs.shape[0], 1)))
        self.rotation_degrees = 0

    def rotate_data(self, degrees):
        self.rotation_degrees += degrees
        self.mapped_rcs = np.dot(self.mapped_rcs, misc.right_rotation_matrix(degrees, degrees=True))
        self.mapped_rcs -= np.tile(self.mapped_rcs.min(axis=0), (self.mapped_rcs.shape[0], 1))
        self.image_shape = map(int, self.mapped_rcs.max(axis=0) + 1)
        return self.image_shape

    def image(self):
        image = np.zeros(self.image_shape)
        image[self.mapped_rcs.astype(np.int)[:,0], self.mapped_rcs.astype(np.int)[:,1]] = 1
        sigma = 0.25 / self.um_per_pixel  # Clusters have stdev ~= 0.25 um
        image = gaussian_filter(image, sigma)
        return image

    def fft_align_with_im(self, image_data):
        im_data_im_shapes = set(a.shape for a in image_data.all_ffts.values())
        assert len(im_data_im_shapes) <= 2, im_data_im_shapes

        # Make the ffts
        fq_image = self.image()
        fq_im_fft_given_shape = {}
        for shape in im_data_im_shapes:
            padded_fq_im = misc.pad_to_size(fq_image, shape)
            fq_im_fft_given_shape[shape] = np.fft.fft2(padded_fq_im)

        # Align
        self.best_max_corr = float('-inf')
        for im_key, im_data_fft in image_data.all_ffts.items():
            fq_im_fft = fq_im_fft_given_shape[im_data_fft.shape]
            cross_corr = abs(np.fft.ifft2(np.conj(fq_im_fft) * im_data_fft))
            max_corr = cross_corr.max()
            max_idx = misc.max_2d_idx(cross_corr)

            if max_corr > self.best_max_corr:
                self.best_im_key = im_key
                self.best_max_corr = max_corr
                self.align_tr = np.array(max_idx) - fq_image.shape
        #print 'Result:', self.key, self.best_im_key, self.best_max_corr, self.align_tr
        return self.key, self.best_im_key, self.best_max_corr, self.align_tr

    def get_new_aligned_rcs(self, new_fq_w=None, new_degree_rot=0, new_tr=(0, 0)):
        """Returns aligned rcs. Only works when image need not be flipped or rotated."""
        if new_fq_w is None:
            new_fq_w = self.w
        aligned_rcs = deepcopy(self.mapped_rcs)
        aligned_rcs = np.dot(aligned_rcs, misc.right_rotation_matrix(new_degree_rot, degrees=True))
        aligned_rcs -= np.tile(aligned_rcs.min(axis=0), (aligned_rcs.shape[0], 1))
        aligned_rcs *= float(new_fq_w) / self.w
        aligned_rcs += np.tile(self.align_tr + new_tr, (aligned_rcs.shape[0], 1))
        return aligned_rcs

    def set_aligned_rcs(self):
        self.aligned_rcs = self.get_new_aligned_rcs()

    def set_aligned_rcs_given_transform(self, lbda, theta, offset):
        """Performs transform calculated in FastqImageCorrelator.least_squares_mapping."""
        A = np.zeros((2*len(self.rcs), 4))
        for i, pt in enumerate(self.rcs):
            xir, yir = pt
            A[2*i, :]   = [xir, -yir, 1, 0]
            A[2*i+1, :] = [yir,  xir, 0, 1]

        x = np.array([lbda * np.cos(theta),
                      lbda * np.sin(theta),
                      offset[0],
                      offset[1]])

        # First update w since it depends on previous scale setting
        self.w = lbda * float(self.w) / self.scale
        #self.w = (self.rcs[:, 0].max() - self.rcs[:, 0].min()) * lbda

        self.scale = lbda
        self.rotation = theta
        self.rotation_degrees = theta * 180.0 / np.pi
        self.offset = offset

        self.aligned_rcs = np.dot(A, x).reshape((len(self.rcs), 2))

    def set_correlation(self, im):
        """Sets alignment correlation. Only works when image need not be flipped or rotated."""
        self.best_max_corr =  sum(im[pt[0], pt[1]] for pt in self.aligned_rcs.astype(np.int)
                                  if 0 <= pt[0] < im.shape[0] and 0 <= pt[1] < im.shape[1])

    def set_snr(self, snr):
        self.snr = snr

    def set_snr_with_control_corr(self, control_corr):
        self.snr = self.best_max_corr / control_corr

    def plot_convex_hull(self, rcs=None, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        if rcs is None:
            rcs = self.aligned_rcs
        hull = ConvexHull(rcs)
        ax.plot(rcs[hull.vertices, 1], rcs[hull.vertices, 0], label=self.key)

    def plot_aligned_rcs(self, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(self.aligned_rcs[:, 1], self.aligned_rcs[:, 0], '.', **kwargs)
