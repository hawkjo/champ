from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.optimize import minimize
import imreg
import misc 


class FastqTileRCs(object):
    """A class for fastq tile coordinates."""
    def __init__(self, key, tile):
        self.key = key
        self.rcs = tile

    def set_fastq_image_data(self, offset, scale, scaled_dims, w, force=False, verbose=True):
        self.offset = offset
        self.scale = scale
        self.image_shape = scaled_dims
        self.w = w  # width in um
        self.mapped_rcs = scale * (self.rcs + np.tile(offset, (self.rcs.shape[0], 1)))
        self.rotation = 0

    def rotate_data(self, degrees):
        self.rotation += degrees
        self.mapped_rcs = np.dot(self.mapped_rcs, misc.right_rotation_matrix(degrees, degrees=True))
        self.mapped_rcs -= np.tile(self.mapped_rcs.min(axis=0), (self.mapped_rcs.shape[0], 1))
        self.image_shape = self.mapped_rcs.max(axis=0) + 1
        return self.image_shape

    def image(self):
        image = np.zeros(self.image_shape)
        image[self.mapped_rcs.astype(np.int)[:,0], self.mapped_rcs.astype(np.int)[:,1]] = 1
        return image

    def imreg_align_with_im(self, im):
        fq_image = self.image()
        edge_len = misc.next_power_of_2(np.r_[fq_image.shape, im.shape].max())
        sq_fq_im = misc.pad_to_size(fq_image, (edge_len, edge_len))

        self.max_score = float('-inf')
        for flip in [False, True]:
            if flip:
                im = np.fliplr(im)
            sq_im = misc.pad_to_size(im, (edge_len, edge_len))
            fq_match_im, scale, rot, tr = imreg.similarity(sq_im, sq_fq_im)
            score = (sq_im * fq_match_im).sum()

            if score > self.max_score:
                self.max_score = score
                self.best_match_im = fq_match_im
                self.align_scale = scale
                self.align_rot = rot
                self.align_tr = tr
        print self.key, score, scale, rot, tr

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

    def correlation(self, im, new_fq_w=None, new_degree_rot=0, new_tr=(0, 0)):
        """Returns alignment correlation. Only works when image need not be flipped or rotated."""
        return sum(im[pt[0], pt[1]] for pt in self.get_new_aligned_rcs(new_fq_w, new_degree_rot, new_tr)
                   if 0 <= pt[0] < im.shape[0] and 0 <= pt[1] < im.shape[1])

    def optimize_alignment(self, im):
        def neg_corr(v):
            return - self.correlation(im, v[0], v[1], v[2:3])
        v0 = [self.w, 0, 0, 0]
        methods = ['Nelder-Mead',
                   'Powell',
                   'CG',
                   'BFGS',
                   'Newton-CG',
                   'Anneal',
                   'L-BFGS-B',
                   'TNC',
                   'COBYLA',
                   'SLSQP',
                   'dogleg',
                   'trust-ncg']
        res = minimize(neg_corr, v0, method=methods[0])
        return res

    def plot_convex_hull(self, rcs=None, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        if rcs is None:
            rcs = self.mapped_rcs
        hull = ConvexHull(rcs)
        ax.plot(rcs[hull.vertices, 1], rcs[hull.vertices, 0], label=self.key)
