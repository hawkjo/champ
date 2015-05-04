"""
A space to create commonly used functions and classes for image processing.
"""
from copy import deepcopy
import itertools
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import convolve2d
from scipy.spatial import ConvexHull
from scipy import ndimage
from pathos.multiprocessing import ProcessingPool
from multiprocessing import Pool
import imreg
import local_config


def signal_hist_and_func(pyplot_hist, plot_title='', plot_curves=True):
    """Determines empirical background noise curve and returns (good signal)/(all signal)
    interpolation function."""

    hy, bin_edges, _ = pyplot_hist
    hx = (bin_edges[:-1] + bin_edges[1:])/2

    mode_y = max(hy[:-1])
    mode_ind = list(hy).index(mode_y)

    gy = np.zeros(hx.shape)
    gy[:mode_ind+1] = hy[:mode_ind+1]
    gy[mode_ind+1 : 2*mode_ind+1] = hy[mode_ind-1::-1]

    sig_y = hy - gy
    sig_y[sig_y < 0] = 0

    ratio = sig_y / hy

    last_zero_ind = np.where(ratio == 0)[0][-1]
    first_one_ind = np.where(ratio >= 1)[0][0]
    ratio[:last_zero_ind] = 0.0
    ratio[first_one_ind:] = 1.0
    l_bnd = hx[last_zero_ind]
    u_bnd = hx[first_one_ind]
    print 'Non-trivial cdf range: %f - %f' % (l_bnd, u_bnd)

    delta_x = np.mean(hx[1:]-hx[:-1])
    print 'delta_x:', delta_x

    num_ext_points = 10
    extended_x = np.r_[[0], hx, 1.1*hx[-1]]
    extended_ratio = np.r_[[0], ratio, [1]]

    ratio_f_interp = interp1d(extended_x, extended_ratio, kind='cubic')

    if plot_curves:
        fig = plt.figure(figsize=(10,10))
        plt.plot(hx, hy, label='Data')
        plt.plot(hx, gy, label='Noise')
        plt.plot(hx, sig_y, label='Signal')
        plt.plot(hx, ratio, label='Ratio')
        plt.plot(extended_x, ratio_f_interp(extended_x), label='Ratio Interp', linewidth=3)
        plt.legend()
        if plot_title:
            plt.title(plot_title)

    return ratio_f_interp


def next_power_of_2(x):
    return 1<<(int(np.ceil(x))-1).bit_length()


def max_2d_idx(a):
    return np.unravel_index(a.argmax(), a.shape)


def pad_to_size(M, size):
    assert len(size) == 2, 'Row and column sizes needed.'
    left_to_pad = size - np.array(M.shape) 
    return np.pad(M, ((0, left_to_pad[0]), (0, left_to_pad[1])), mode='constant')


def right_rotation_matrix(angle, degrees=True):
    if degrees:
        angle *= np.pi/180.0
    sina = np.sin(angle)
    cosa = np.cos(angle)
    return np.array([[cosa, sina],
                     [-sina, cosa]])


class FastqTileXYs(object):
    def __init__(self, key, tile):
        self.key = key
        self.xys = tile

    def set_fastq_image_data(self, offset, scale, scaled_dims, w, force=False, verbose=True):
        self.offset = offset
        self.scale = scale
        self.image_shape = scaled_dims
        self.w = w  # width in um
        self.mapped_xys = scale * (self.xys + np.tile(offset, (self.xys.shape[0], 1)))
        self.rotation = 0

    def rotate_data(self, degrees):
        self.rotation += degrees
        self.mapped_xys = np.dot(self.mapped_xys, right_rotation_matrix(degrees, degrees=True))
        self.mapped_xys -= np.tile(self.mapped_xys.min(axis=0), (self.mapped_xys.shape[0], 1))
        self.image_shape = self.mapped_xys.max(axis=0) + 1
        return self.image_shape

    def image(self):
        image = np.zeros(self.image_shape)
        image[self.mapped_xys.astype(np.int)[:,0], self.mapped_xys.astype(np.int)[:,1]] = 1
        return image

    def imreg_align_with_im(self, im):
        fq_image = self.image()
        edge_len = next_power_of_2(np.r_[fq_image.shape, im.shape].max())
        sq_fq_im = pad_to_size(fq_image, (edge_len, edge_len))

        self.max_score = float('-inf')
        for flip in [False, True]:
            if flip:
                im = np.fliplr(im)
            sq_im = pad_to_size(im, (edge_len, edge_len))
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
            padded_fq_im = pad_to_size(fq_image, shape)
            fq_im_fft_given_shape[shape] = np.fft.fft2(padded_fq_im)

        # Align
        self.best_max_corr = float('-inf')
        for im_key, im_data_fft in image_data.all_ffts.items():
            fq_im_fft = fq_im_fft_given_shape[im_data_fft.shape]
            cross_corr = abs(np.fft.ifft2(np.conj(fq_im_fft) * im_data_fft))
            max_corr = cross_corr.max()
            max_idx = max_2d_idx(cross_corr)

            if max_corr > self.best_max_corr:
                self.best_im_key = im_key
                self.best_max_corr = max_corr
                self.align_tr = np.array(max_idx) - fq_image.shape
        #print 'Result:', self.key, self.best_im_key, self.best_max_corr, self.align_tr
        return self.key, self.best_im_key, self.best_max_corr, self.align_tr

    def aligned_xys(self, new_fq_w=None, new_degree_rot=0, new_tr=(0, 0)):
        """Returns aligned xys. Only works when image need not be flipped or rotated."""
        aligned_xys = deepcopy(self.mapped_xys)
        if new_fq_w:
            aligned_xys *= float(new_fq_w) / self.w
        aligned_xys += np.tile(self.align_tr + new_tr, (aligned_xys.shape[0], 1))
        aligned_xys = np.dot(aligned_xys, right_rotation_matrix(new_degree_rot, degrees=True))
        # Fix rotation and flip offsets and align:
        return aligned_xys

    def correlation(self, im, new_fq_w=None, new_degree_rot=0, new_tr=(0, 0)):
        """Returns alignment correlation. Only works when image need not be flipped or rotated."""
        return sum(im[pt[0], pt[1]] for pt in self.aligned_xys(new_fq_w, new_degree_rot, new_tr)
                   if 0 <= pt[0] < im.shape[0] and 0 <= pt[1] < im.shape[1])

    def plot_convex_hull(self, xys=None, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        if xys is None:
            xys = self.mapped_xys
        hull = ConvexHull(xys)
        ax.plot(xys[hull.vertices, 1], xys[hull.vertices, 0])


class ImageData(object):
    def __init__(self, im, objective):
        assert isinstance(im, np.ndarray), 'Image must be numpy ndarray'
        assert objective in set([40, 60]), 'Accepted objectives are 40 and 60'

        self.im = im
        self.objective = objective
        self.um_per_pixel = 16.0 / self.objective
        self.um_dims = self.um_per_pixel * np.array(self.im.shape)

    def D4_im_given_idx(self, idx):
        flip, rot = idx
        if flip:
            flip_im = np.fliplr(self.im)
        else:
            flip_im = self.im
        return np.rot90(flip_im, k=(rot/90))

    def iterate_D4_images(self):
        for idx in self.iterate_D4_idxs():
            yield idx, self.D4_im_given_idx(idx)

    def iterate_D4_idxs(self):
        for flip in [0, 1]:
            for rot in [0, 90, 180, 270]:
                yield (flip, rot)

    def D4_ffts(self, padding=(0, 0), processors=1, force=False):
        """Makes images and ffts of all flips and 90 degree rotations (i.e. D4)"""
        if hasattr(self, 'all_ffts') and self.all_ffts and not force:
            return
        self.fft_padding = padding
        pool = ProcessingPool(processors)
        ffts = pool.map(self.single_fft, self.iterate_D4_idxs())
        self.all_ffts = {idx: fft for idx, fft in zip(self.iterate_D4_idxs(), ffts)}
        del pool

    def set_single_fft(self, idx, padding):
        self.fft_padding = padding
        self.all_ffts = {idx: self.single_fft(idx)}

    def single_fft(self, idx):
        im = self.D4_im_given_idx(idx)
        totalx, totaly = np.array(self.fft_padding)+np.array(im.shape)
        w = next_power_of_2(totalx)
        h = next_power_of_2(totaly)
        padded_im = np.pad(im,
                           ((self.fft_padding[0], w-totalx), (self.fft_padding[1], h-totaly)),
                           mode='constant')
        return np.fft.fft2(padded_im)


class FastqImageCorrelator(object):
    def __init__(self, project_name):
        self.project_name = project_name
        self.fastq_tiles = {}
        self.fastq_tiles_list = []
        self.fastq_tiles_keys = []
        self.image_data = None
        self.w_fq_tile_min = 895  # um
        self.w_fq_tile_max = 937  # um
        self.w_fq_tile = 937  # um

    def load_phiX(self):
        for key, tile in local_config.phiX_xys_given_project_name(self.project_name).items():
            self.fastq_tiles[key] = FastqTileXYs(key, tile)
        self.fastq_tiles_keys = [key for key, tile in sorted(self.fastq_tiles.items())]
        self.fastq_tiles_list = [tile for key, tile in sorted(self.fastq_tiles.items())]

    def set_image_data(self, im):
        assert isinstance(im, ImageData), 'Object passed to set_image_data must be ImageData object.'
        self.image_data = im

    def set_image_data_from_ndarray(self, im, objective):
        self.image_data = ImageData(im, objective)

    def load_image_data_from_fpath(self, fpath, objective):
        self.data_im = ImageData(np.load(fpath), objective)

    def set_fastq_tile_mappings(self):
        """Calculate parameters for mapping fastq tiles for ffts."""
        assert self.image_data is not None, 'No image data loaded.'
        assert self.fastq_tiles != {}, 'No fastq data loaded.'

        self.all_data = np.concatenate([tile.xys for key, tile in self.fastq_tiles.items()])
    
        x_min, y_min = self.all_data.min(axis=0)
        x_max, y_max = self.all_data.max(axis=0)
    
        self.fq_im_offset = np.array([-x_min, -y_min])
        self.fq_im_scale = (float(self.w_fq_tile) / (x_max-x_min)) / self.image_data.um_per_pixel
        self.fq_im_scaled_maxes = self.fq_im_scale * np.array([x_max-x_min, y_max-y_min])
        self.fq_im_scaled_dims = (self.fq_im_scaled_maxes + [1, 1]).astype(np.int)

    def set_all_fastq_image_data(self, verbose=True):
        for key, tile in self.fastq_tiles.items():
            tile.set_fastq_image_data(self.fq_im_offset,
                                      self.fq_im_scale,
                                      self.fq_im_scaled_dims,
                                      self.w_fq_tile,
                                      verbose=verbose)

    def rotate_all_fastq_data(self, degrees):
        im_shapes = []
        for tile in self.fastq_tiles_list:
            im_shapes.append(tile.rotate_data(degrees))
        self.fq_im_scaled_dims = np.array(im_shapes).max(axis=0)
        for tile in self.fastq_tiles_list:
            tile.image_shape = self.fq_im_scaled_dims


    def imreg_align(self):
        for key, tile in sorted(self.fastq_tiles.items()):
            tile.imreg_align_with_im(self.image_data.im)

    def fft_align_tile(self, tile):
        return tile.fft_align_with_im(self.image_data)

    def fft_align_tile_with_im(self, tile):
        im_data_im_shapes = set(a.shape for a in self.image_data.all_ffts.values())
        assert len(im_data_im_shapes) <= 2, im_data_im_shapes

        # Make the ffts
        fq_image = tile.image()
        fq_im_fft_given_shape = {}
        for shape in im_data_im_shapes:
            padded_fq_im = pad_to_size(fq_image, shape)
            fq_im_fft_given_shape[shape] = np.fft.fft2(padded_fq_im)

        # Align
        best_max_corr = float('-inf')
        for im_key, im_data_fft in self.image_data.all_ffts.items():
            fq_im_fft = fq_im_fft_given_shape[im_data_fft.shape]
            cross_corr = abs(np.fft.ifft2(np.conj(fq_im_fft) * im_data_fft))
            max_corr = cross_corr.max()
            max_idx = max_2d_idx(cross_corr)

            if max_corr > best_max_corr:
                best_im_key = im_key
                best_max_corr = max_corr
                align_tr = np.array(max_idx) - fq_image.shape
        #print 'Result:', tile.key, best_im_key, best_max_corr, align_tr
        return tile.key, best_im_key, best_max_corr, align_tr

    def fft_align(self, processors, recalc_fft=True, verbose=True):
        if verbose:
            print 'Set fastq tile mappings'
        self.set_fastq_tile_mappings()
        if verbose:
            print 'Image D4 ffts'
        self.image_data.D4_ffts(padding=self.fq_im_scaled_dims,
                                processors=processors,
                                force=recalc_fft)
        if verbose:
            print 'Fastq images and ffts'
        self.set_all_fastq_image_data(verbose=True)
        if verbose:
            print 'Aligning'
        self.best_corr = 0
        self.max_corrs = []
        for tile in self.fastq_tiles_list:
            fq_key, im_key, max_corr, align_tr = self.fft_align_tile(tile)
            self.max_corrs.append(max_corr)
            if max_corr > self.best_corr:
                self.best_corr = max_corr
                self.best_fq_key = fq_key
                self.best_im_key = im_key
                self.best_align_tr = align_tr
        if verbose:
            print 'Best result:', self.best_corr, self.best_fq_key, self.best_im_key, self.best_align_tr

    def show_alignment(self, fq_key, im_key, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        im = self.image_data.D4_im_given_idx(im_key)
        ax.imshow(im, cmap=plt.get_cmap('Blues'), label='Ilya')
        v = plt.axis()
        fq_tile = self.fastq_tiles[fq_key]
        fq_im = fq_tile.image()[-fq_tile.align_tr[0]:, -fq_tile.align_tr[1]:]
        fq_im = ndimage.filters.gaussian_filter(fq_im, 3)
        fq_im[fq_im < 0.01] = None
        plt.imshow(fq_im, cmap=plt.get_cmap('Reds'), alpha=0.5, label='Fastq')
        #aligned_xys = self.fastq_tiles[fq_key].aligned_xys()
        #ax.plot(aligned_xys[:, 0], aligned_xys[:, 1], 'r.', alpha=0.5)
        plt.axis(v)
