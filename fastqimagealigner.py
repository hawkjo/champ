from copy import deepcopy
from fastqtilercs import FastqTileRCs
from imagedata import ImageData
from itertools import izip
import logging
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from misc import pad_to_size, max_2d_idx, AlignmentStats
import numpy as np
import random
import reads
from scipy import ndimage
import scipy.optimize
from scipy.spatial import KDTree
import sextraction
from sklearn.mixture import GMM
import time

log = logging.getLogger(__name__)


class FastqImageAligner(object):
    """A class to find the alignment of fastq data and image data."""
    def __init__(self, chip_id, file_structure):
        self.chip_id = chip_id
        self.file_structure = file_structure
        self.fastq_tiles = {}
        self.fastq_tiles_list = []
        self.fastq_tiles_keys = []
        self.image_data = ImageData()
        self.w_fq_tile_min = 895  # um
        self.w_fq_tile_max = 937  # um
        self.fq_w = 927  # um

    def load_phiX(self):
        tile_data = reads.phix_read_names(self.chip_id, self.file_structure)
        self.load_reads(tile_data=tile_data)

    def load_all_reads(self, tile_keys=None):
        tile_data = reads.all_read_names(self.chip_id, self.file_structure)
        if tile_keys is None:
            self.load_reads(tile_data)
        else:
            self.load_reads({tile_key: read_names for tile_key, read_names in tile_data.items() if tile_key in tile_keys})

    def load_reads(self, tile_data):
        for tile_key, read_names in tile_data.items():
            self.fastq_tiles[tile_key] = FastqTileRCs(tile_key, read_names)
        self.fastq_tiles_keys = [tile_key for tile_key, tile in sorted(self.fastq_tiles.items())]
        self.fastq_tiles_list = [tile for tile_key, tile in sorted(self.fastq_tiles.items())]

    def all_reads_fic_from_aligned_fic(self, other_fic, tile_data=None):
        if tile_data is None:
            self.load_all_reads(tile_keys=[tile.key for tile in other_fic.hitting_tiles])
        else:
            self.load_reads(tile_data)
        self.image_data = deepcopy(other_fic.image_data)
        self.fq_w = other_fic.fq_w
        self.set_fastq_tile_mappings()
        self.set_all_fastq_image_data()
        self.hitting_tiles = [self.fastq_tiles[tile.key] for tile in other_fic.hitting_tiles]
        self.sexcat = other_fic.sexcat

        for other_tile in other_fic.hitting_tiles:
            tile = self.fastq_tiles[other_tile.key]
            tile.set_aligned_rcs_given_transform(other_tile.scale,
                                                 other_tile.rotation,
                                                 other_tile.offset)

    def fic_given_alignment(self, tile_data, tile_key, scale, fq_w, rotation, rc_offset):
        self.load_reads(tile_data)
        self.hitting_tiles = []
        self.set_tile_alignment(tile_key, scale, fq_w, rotation, rc_offset)

    def set_tile_alignment(self, tile_key, scale, fq_w, rotation, rc_offset):
        if self.fastq_tiles[tile_key] not in self.hitting_tiles:
            self.hitting_tiles.append(self.fastq_tiles[tile_key])
        self.fq_w = fq_w
        self.set_fastq_tile_mappings()
        self.set_all_fastq_image_data()
        tile = self.fastq_tiles[tile_key]
        tile.set_aligned_rcs_given_transform(scale, rotation, rc_offset)

    def alignment_from_alignment_file(self, fpath):
        self.hitting_tiles = []
        astats = AlignmentStats(fpath)
        for i in range(astats.numtiles):
            self.set_tile_alignment(astats.tile[i],
                                    astats.scaling[i],
                                    astats.tile_width[i],
                                    astats.rotation[i] * np.pi / 180,
                                    astats.rc_offset[i]
                                   )
        
    def set_image_data(self, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], ImageData):
            self.image_data = args[0]
        else:
            self.image_data = ImageData(*args, **kwargs)

    def set_fastq_tile_mappings(self):
        """Calculate parameters for mapping fastq tiles for ffts."""
        assert self.image_data is not None, 'No image data loaded.'
        assert self.fastq_tiles != {}, 'No fastq data loaded.'

        self.all_data = np.concatenate([tile.rcs for tile in self.fastq_tiles.values()])
        print("all data len", len(self.all_data))
        x_min, y_min = self.all_data.min(axis=0)
        x_max, y_max = self.all_data.max(axis=0)
        print("xmin, ymin", x_min, y_min)
        print("x_max, y_max ", x_max, y_max)
        self.fq_im_offset = np.array([-x_min, -y_min])
        self.fq_im_scale = (float(self.fq_w) / (x_max-x_min)) / self.image_data.um_per_pixel
        self.fq_im_scaled_maxes = self.fq_im_scale * np.array([x_max-x_min, y_max-y_min])
        self.fq_im_scaled_dims = (self.fq_im_scaled_maxes + [1, 1]).astype(np.int)
        print("the four in set fastq tile mappings", self.fq_im_offset, self.fq_im_scale, self.fq_im_scaled_maxes, self.fq_im_scaled_dims)

    def set_all_fastq_image_data(self):
        print("tile shapes in set all fastq image data")
        for tile in self.fastq_tiles.values():
            tile.set_fastq_image_data(self.fq_im_offset,
                                      self.fq_im_scale,
                                      self.fq_im_scaled_dims,
                                      self.fq_w)
            print(tile.image_shape)

    def rotate_all_fastq_data(self, degrees):
        print("rotating fastq data")
        im_shapes = [tile.rotate_data(degrees) for tile in self.fastq_tiles_list]
        self.fq_im_scaled_dims = np.array(im_shapes).max(axis=0)
        for tile in self.fastq_tiles_list:
            tile.image_shape = self.fq_im_scaled_dims
            print(tile.image_shape)

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
        best_im_key = None
        align_tr = None
        for im_key, im_data_fft in self.image_data.all_ffts.items():
            fq_im_fft = fq_im_fft_given_shape[im_data_fft.shape]
            cross_corr = abs(np.fft.ifft2(np.conj(fq_im_fft) * im_data_fft))
            max_corr = cross_corr.max()
            max_idx = max_2d_idx(cross_corr)

            if max_corr > best_max_corr:
                best_im_key = im_key
                best_max_corr = max_corr
                align_tr = np.array(max_idx) - fq_image.shape
        if best_im_key is None:
            raise ValueError("Unable to align tile")
        return tile.key, best_im_key, best_max_corr, align_tr

    def fft_align(self, processors, recalc_fft=True):
        log.debug('Set fastq tile mappings')
        self.set_fastq_tile_mappings()
        log.debug('Image D4 ffts')
        self.image_data.D4_ffts(padding=self.fq_im_scaled_dims,
                                processors=processors,
                                force=recalc_fft)
        log.debug('Fastq images and ffts')
        self.set_all_fastq_image_data()
        log.debug('Aligning')
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
        log.debug('Best result: %s, %s, %s, %s' % (self.best_corr, self.best_fq_key, self.best_im_key, self.best_align_tr))

    def show_alignment(self, fq_key, im_key, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        im = self.image_data.D4_im_given_idx(im_key)
        ax.matshow(im, cmap=plt.get_cmap('Blues'), label='Ilya')
        v = plt.axis()
        fq_tile = self.fastq_tiles[fq_key]
        fq_im = fq_tile.image()[-fq_tile.align_tr[0]:, -fq_tile.align_tr[1]:]
        fq_im = ndimage.filters.gaussian_filter(fq_im, 3)
        fq_im[fq_im < 0.01] = None
        plt.matshow(fq_im, cmap=plt.get_cmap('Reds'), alpha=0.5, label='Fastq')
        plt.axis(v)

    def set_sexcat(self, sexcat):
        assert isinstance(sexcat, sextraction.Sextraction)
        self.sexcat = sexcat

    def set_sexcat_from_file(self, fpath):
        self.sexcat = sextraction.Sextraction(fpath)

    def find_hitting_tiles(self, possible_tile_keys, snr_thresh=1.2):
        possible_tiles = [self.fastq_tiles[key] for key in possible_tile_keys]
        impossible_tiles = [tile for tile in self.fastq_tiles.values() if tile not in possible_tiles]
        control_tiles = random.sample(impossible_tiles, 3)

        self.image_data.set_single_fft((0, 0), padding=self.fq_im_scaled_dims)
        self.control_corr = 0
        for control_tile in control_tiles:
            corr = self.fft_align_tile(control_tile)[2]
            if corr > self.control_corr:
                self.control_corr = corr

        self.hitting_tiles = []
        for tile in possible_tiles:
            tile_key, _, max_corr, align_tr = self.fft_align_tile(tile)
            if max_corr > snr_thresh * self.control_corr:
                tile.set_aligned_rcs()
                tile.set_snr(max_corr / self.control_corr)
                self.hitting_tiles.append(tile)

    def find_points_in_frame(self, consider_tiles='all'):
        self.aligned_rcs_in_frame = []
        self.rcs_in_frame = []
        im_shape = self.image_data.im.shape

        if consider_tiles == 'all':
            considered_tiles = self.hitting_tiles
        elif consider_tiles == 'best':
            best_corr = 0
            for tile in self.hitting_tiles:
                if tile.best_max_corr > best_corr:
                    considered_tiles = [tile]
                    best_corr = tile.best_max_corr
        elif isinstance(consider_tiles, FastqTileRCs):
            considered_tiles = [consider_tiles]
        elif isinstance(consider_tiles, int):
            considered_tiles = [self.hitting_tiles[consider_tiles]]
        elif isinstance(consider_tiles, list):
            if np.all([isinstance(el, FastqTileRCs) for el in consider_tiles]):
                considered_tiles = consider_tiles
            elif np.all([isinstance(el, int) for el in consider_tiles]):
                considered_tiles = [self.hitting_tiles[i] for i in consider_tiles]
            else:
                raise ValueError('A consider_tiles list must either be all tiles or all indices')
        else:
            raise ValueError('Invalid consider_tiles parameter.')

        for tile in considered_tiles:
            rcs = tile.rcs.astype(np.int)
            for i, pt in enumerate(tile.aligned_rcs):
                if 0 <= pt[0] < im_shape[0] and 0 <= pt[1] < im_shape[1]:
                    self.aligned_rcs_in_frame.append(pt)
                    self.rcs_in_frame.append((tile.key, rcs[i]))
        self.aligned_rcs_in_frame = np.array(self.aligned_rcs_in_frame)

    def hit_dists(self, hits):
        return [self.single_hit_dist(hit) for hit in hits]

    def single_hit_dist(self, hit):
        return np.linalg.norm(self.sexcat.point_rcs[hit[0]] - self.aligned_rcs_in_frame[hit[1]])

    def remove_longest_hits(self, hits, pct_thresh):
        dists = self.hit_dists(hits)
        thresh = np.percentile(dists, pct_thresh * 100)
        return [hit for hit in hits if self.single_hit_dist(hit) <= thresh]

    def find_hits(self, second_neighbor_thresh=None, consider_tiles='all'):
        # --------------------------------------------------------------------------------
        # Find nearest neighbors
        # --------------------------------------------------------------------------------
        self.find_points_in_frame(consider_tiles)
        sexcat_tree = KDTree(self.sexcat.point_rcs)
        aligned_tree = KDTree(self.aligned_rcs_in_frame)

        # All indices are in the order (sexcat_idx, aligned_in_frame_idx)
        sexcat_to_aligned_idxs = set()
        for i, pt in enumerate(self.sexcat.point_rcs):
            dist, idx = aligned_tree.query(pt)
            sexcat_to_aligned_idxs.add((i, idx))

        aligned_to_sexcat_idxs_rev = set()
        for i, pt in enumerate(self.aligned_rcs_in_frame):
            dist, idx = sexcat_tree.query(pt)
            aligned_to_sexcat_idxs_rev.add((idx, i))

        # --------------------------------------------------------------------------------
        # Find categories of hits
        # --------------------------------------------------------------------------------
        mutual_hits = sexcat_to_aligned_idxs & aligned_to_sexcat_idxs_rev
        non_mutual_hits = sexcat_to_aligned_idxs ^ aligned_to_sexcat_idxs_rev

        sexcat_in_non_mutual = set(i for i, j in non_mutual_hits)
        aligned_in_non_mutual = set(j for i, j in non_mutual_hits)
        exclusive_hits = set((i, j) for i, j in mutual_hits if i not in
                             sexcat_in_non_mutual and j not in aligned_in_non_mutual)

        # --------------------------------------------------------------------------------
        # Recover good non-exclusive mutual hits. 
        # --------------------------------------------------------------------------------
        # If the distance to second neighbor is too close, that suggests a bad peak call combining
        # two peaks into one. Filter those out with a gaussian-mixture-model-determined threshold.
        if self.image_data.objective == 60:
            # Value decided by observation of our data. May vary with equipment.
            good_hit_thresh = 5
        else:
            good_hit_thresh = np.percentile(self.hit_dists(exclusive_hits), 95)

        if second_neighbor_thresh is not None:
            assert second_neighbor_thresh > 0
            self.second_neighbor_thresh = second_neighbor_thresh
        else:
            self.second_neighbor_thresh = 2 * good_hit_thresh

        exclusive_hits = set(hit for hit in exclusive_hits
                             if self.single_hit_dist(hit) <= good_hit_thresh)

        good_mutual_hits = set()
        for i, j in (mutual_hits - exclusive_hits):
            if self.hit_dists([(i, j)])[0] > good_hit_thresh:
                continue
            third_wheels = [tup for tup in non_mutual_hits if i == tup[0] or j == tup[1]]
            if min(self.hit_dists(third_wheels)) > self.second_neighbor_thresh:
                good_mutual_hits.add((i, j))
        bad_mutual_hits = mutual_hits - exclusive_hits - good_mutual_hits

        # --------------------------------------------------------------------------------
        # Test that the four groups form a partition of all hits and finalize
        # --------------------------------------------------------------------------------
        assert (non_mutual_hits | bad_mutual_hits | good_mutual_hits | exclusive_hits
                == sexcat_to_aligned_idxs | aligned_to_sexcat_idxs_rev
                and len(non_mutual_hits) + len(bad_mutual_hits)
                + len(good_mutual_hits) + len(exclusive_hits)
                == len(sexcat_to_aligned_idxs | aligned_to_sexcat_idxs_rev))

        self.non_mutual_hits = non_mutual_hits
        self.mutual_hits = mutual_hits
        self.bad_mutual_hits = bad_mutual_hits
        self.good_mutual_hits = good_mutual_hits
        self.exclusive_hits = exclusive_hits

        log.info('Non-mutual hits: %s' % len(non_mutual_hits))
        log.info('Mutual hits: %s' % len(mutual_hits))
        log.info('Bad mutual hits: %s' % len(bad_mutual_hits))
        log.info('Good mutual hits: %s' % len(good_mutual_hits))
        log.info('Exclusive hits: %s' % len(exclusive_hits))

    def least_squares_mapping(self, hit_type='exclusive', pct_thresh=0.9, min_hits=50):
        """least_squares_mapping(self, hit_type='exclusive')

        "Input": set of tuples of (sexcat_idx, in_frame_idx) mappings.

        "Output": scaling lambda, rotation theta, x_offset, y_offset, and aligned_rcs

        We here solve the matrix least squares equation Ax = b, where

                [ x0r -y0r 1 0 ]
                [ y0r  x0r 0 1 ]
            A = [ x1r -y1r 1 0 ]
                [ y1r  x1r 0 1 ]
                      . . .
                [ xnr -ynr 1 0 ]
                [ ynr  xnr 0 1 ]

        and

            b = [ x0s y0s x1s y1s . . . xns yns ]^T

        The r and s subscripts indicate rcs and sexcat coords.

        The interpretation of x is then given by

            x = [ alpha beta x_offset y_offset ]^T

        where
            alpha = lambda cos(theta), and
            beta = lambda sin(theta)

        This system of equations is then finally solved for lambda and theta.
        """
        def get_hits(hit_type):
            if isinstance(hit_type, str):
                hit_type = [hit_type]
            hits = []
            for ht in hit_type:
                hits.extend(getattr(self, ht + '_hits'))
            return hits 

        for tile in self.hitting_tiles:
            self.find_hits(consider_tiles=tile)

            # Reminder: All indices are in the order (sexcat_idx, in_frame_idx)
            hits = self.remove_longest_hits(get_hits(hit_type), pct_thresh)
            assert len(hits) > min_hits, 'Too few hits for least squares mapping: {0}'.format(len(hits))
            A = np.zeros((2 * len(hits), 4))
            b = np.zeros((2 * len(hits),))
            for i, (sexcat_idx, in_frame_idx) in enumerate(hits):
                tile_key, (xir, yir) = self.rcs_in_frame[in_frame_idx]
                A[2*i, :] = [xir, -yir, 1, 0]
                A[2*i+1, :] = [yir,  xir, 0, 1]

                xis, yis = self.sexcat.point_rcs[sexcat_idx]
                b[2*i] = xis
                b[2*i+1] = yis

            alpha, beta, x_offset, y_offset = np.linalg.lstsq(A, b)[0]
            offset = np.array([x_offset, y_offset])
            theta = np.arctan2(beta, alpha)
            lbda = alpha / np.cos(theta)

            tile.set_aligned_rcs_given_transform(lbda, theta, offset)
            tile.set_correlation(self.image_data.im)
            if hasattr(self, 'control_corr'):
                tile.set_snr_with_control_corr(self.control_corr)

    def align(self, possible_tile_keys, rotation_est, fq_w_est=927, snr_thresh=1.2,
              hit_type=('exclusive', 'good_mutual'), min_hits=15):

        start_time = time.time()
        self.fq_w = fq_w_est
        self.set_fastq_tile_mappings()
        self.set_all_fastq_image_data()
        self.rotate_all_fastq_data(rotation_est)
        log.debug('Prep time: %.3f seconds' % (time.time() - start_time))

        start_time = time.time()
        self.find_hitting_tiles(possible_tile_keys, snr_thresh)
        log.debug('Rough alignment time: %.3f seconds' % (time.time() - start_time))

        start_time = time.time()
        if not self.hitting_tiles:
            raise RuntimeError('Alignment not found')
        self.least_squares_mapping(hit_type, min_hits=min_hits)
        log.debug('Precision alignment time: %.3f seconds' % (time.time() - start_time))

        start_time = time.time()
        self.find_hits()
        log.debug('Hit finding time: %.3f seconds' % (time.time() - start_time))
        
    def precision_align_only(self, hit_type=('exclusive', 'good_mutual'), min_hits=15):
        start_time = time.time()
        if not self.hitting_tiles:
            raise RuntimeError('Alignment not found')
        self.least_squares_mapping(hit_type, min_hits=min_hits)
        log.debug('Precision alignment time: %.3f seconds' % (time.time() - start_time))
        
        start_time = time.time()
        self.find_hits()
        log.debug('Hit finding time: %.3f seconds' % (time.time() - start_time))
        
    def gmm_thresh(self, dists):
        self.gmm = GMM(2)
        self.gmm.fit(dists)
        lower_idx = self.gmm.means_.argmin()
        higher_idx = 1 - lower_idx
        lower_mean = self.gmm.means_[lower_idx]
        good_posterior_thresh_pct = 0.99
        f = lambda x: self.gmm.predict_proba([x])[0][higher_idx] - good_posterior_thresh_pct
        self.second_neighbor_thresh = scipy.optimize.brentq(f, lower_mean, max(dists))

    def plot_hit_hists(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 8))
        non_mut_dists = self.hit_dists(self.non_mutual_hits)
        bins = np.linspace(0, max(non_mut_dists), 50)

        if non_mut_dists:
            ax.hist(non_mut_dists, bins, label='Non-mutual hits', normed=True, histtype='step')
        if self.bad_mutual_hits:
            ax.hist(self.hit_dists(self.bad_mutual_hits), bins, label='Bad mutual hits', normed=True, histtype='step')
        if self.good_mutual_hits:
            ax.hist(self.hit_dists(self.good_mutual_hits), bins, label='Good mutual hits', normed=True, histtype='step')
        if self.exclusive_hits:
            ax.hist(self.hit_dists(self.exclusive_hits), bins, label='Exclusive hits', normed=True, histtype='step')
        ax.legend()
        ax.set_title('%s Nearest Neighbor Distance Distributions' % self.image_data.bname)
        return ax

    def plot_threshold_gmm(self, axs=None, force=False):
        if axs is None:
            fig, axs = plt.subplots(1, 2, figsize=(15, 6))
        non_mut_dists = self.hit_dists(self.non_mutual_hits)
        if not hasattr(self, 'gmm') and force:
            self.gmm_thresh(non_mut_dists)
        xs = np.linspace(0, max(non_mut_dists), 200)
        posteriors = self.gmm.predict_proba(xs)
        pdf = np.exp(self.gmm.score_samples(xs)[0])

        axs[0].hist(non_mut_dists, 40, histtype='step', normed=True, label='Data')
        axs[0].plot(xs, pdf, label='PDF')
        ylim = axs[0].get_ylim()
        axs[0].plot([self.second_neighbor_thresh, self.second_neighbor_thresh], ylim, 'g--', label='Threshold')
        axs[0].set_title('%s GMM PDF of Non-mutual hits' % self.image_data.bname)
        axs[0].legend()
        axs[0].set_ylim(ylim)

        axs[1].hist(non_mut_dists, 40, histtype='step', normed=True, label='Data')
        axs[1].plot(xs, posteriors, label='Posterior')
        axs[1].plot([self.second_neighbor_thresh, self.second_neighbor_thresh], [0, 1], 'g--', label='Threshold')
        axs[1].set_title('%s GMM Posterior Probabilities' % self.image_data.bname)
        axs[1].legend()
        return axs

    def plot_hits(self, hits, color, ax, kwargs):
        for i, j in hits:
            ax.plot([self.sexcat.point_rcs[i, 1], self.aligned_rcs_in_frame[j, 1]],
                    [self.sexcat.point_rcs[i, 0], self.aligned_rcs_in_frame[j, 0]],
                    color=color, **kwargs)
        return ax

    def plot_all_hits(self, ax=None, im_kwargs=None, line_kwargs=None, fqpt_kwargs=None, sext_kwargs=None,
                     title_kwargs=None, legend_kwargs=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(15, 15))

        kwargs = {'cmap': plt.get_cmap('Blues')}
        kwargs.update(im_kwargs or {})
        ax.matshow(self.image_data.im, **kwargs)

        kwargs = {'color': 'k', 'alpha': 0.3, 'linestyle': '', 'marker': 'o', 'markersize': 3}
        kwargs.update(fqpt_kwargs or {})
        ax.plot(self.aligned_rcs_in_frame[:, 1], self.aligned_rcs_in_frame[:, 0], **kwargs)

        kwargs = {'alpha': 0.6, 'color': 'darkgoldenrod'}
        kwargs.update(sext_kwargs or {})
        self.sexcat.plot_ellipses(ax=ax, **kwargs)

        line_kwargs = line_kwargs or {}
        title_kwargs = title_kwargs or {}
        legend_kwargs = legend_kwargs or {}
        self.plot_hits(self.non_mutual_hits, 'grey', ax, line_kwargs)
        self.plot_hits(self.bad_mutual_hits, 'b', ax, line_kwargs)
        self.plot_hits(self.good_mutual_hits, 'magenta', ax, line_kwargs)
        self.plot_hits(self.exclusive_hits, 'r', ax, line_kwargs)
        ax.set_title('All Hits: %s vs. %s %s\nRot: %s deg, Fq width: %s um, Scale: %s px/fqu, Corr: %s, SNR: %s'
                % (self.image_data.bname,
                   self.chip_id,
                   ','.join(tile.key for tile in self.hitting_tiles),
                   ','.join('%.2f' % tile.rotation_degrees for tile in self.hitting_tiles),
                   ','.join('%.2f' % tile.w for tile in self.hitting_tiles),
                   ','.join('%.5f' % tile.scale for tile in self.hitting_tiles),
                   ','.join('%.1f' % tile.best_max_corr for tile in self.hitting_tiles),
                   ','.join('%.2f' % tile.snr if hasattr(tile, 'snr') else '-' for tile in self.hitting_tiles),
                   ), **title_kwargs)
        ax.set_xlim([0, self.image_data.im.shape[1]])
        ax.set_ylim([self.image_data.im.shape[0], 0])

        grey_line = Line2D([], [], color='grey', label='Non-mutual hits: %d' % (len(self.non_mutual_hits)))
        blue_line = Line2D([], [], color='blue', label='Bad mutual hits: %d' % (len(self.bad_mutual_hits)))
        magenta_line = Line2D([], [], color='magenta', label='Good mutual hits: %d' % (len(self.good_mutual_hits)))
        red_line = Line2D([], [], color='red', label='Exclusive hits: %d' % (len(self.exclusive_hits)))
        sexcat_line = Line2D([], [], color='darkgoldenrod', alpha=0.6, marker='o', markersize=10,
                             label='Sextractor Ellipses: %d' % (len(self.sexcat.point_rcs)))
        fastq_line = Line2D([], [], color='k', alpha=0.3, marker='o', markersize=10,
                            label='Fastq Points: %d' % (len(self.aligned_rcs_in_frame)))
        handles = [grey_line, blue_line, magenta_line, red_line, sexcat_line, fastq_line]
        legend = ax.legend(handles=handles, **legend_kwargs)
        legend.get_frame().set_color('white')
        return ax

    def plot_hit_vectors(self, hit_types=('exclusive',), ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(15, 15))
        colors = {'exclusive': 'r',
                  'good_mutual': 'magenta',
                  'bad_mutual': 'b',
                  'non_mutual': 'grey'}
        for hit_type in hit_types:
            hits = getattr(self, hit_type + '_hits')
            pts = np.array([self.sexcat.point_rcs[i] - self.aligned_rcs_in_frame[j] for i, j in hits])
            ax.plot(pts[:, 1], pts[:, 0], '.', color=colors[hit_type])
        ax.plot([0], [0], 'k*')
        ax.set_aspect(1)
        ylim = ax.get_ylim()
        ax.set_ylim((ylim[1], ylim[0]))
        ax.set_title('{0} {1} Hit Diffs'.format(self.image_data.bname, hit_type.capitalize()))
        ax.set_xlabel('c')
        ax.set_ylabel('r')
        return ax

    def output_intensity_results(self, out_fpath):
        hit_given_aligned_idx = {}
        for hit_type in ('non_mutual', 'bad_mutual', 'good_mutual', 'exclusive'):
            for i, j in getattr(self, hit_type + '_hits'):
                hit_given_aligned_idx[j] = (hit_type, (i, j))

        hit_given_rcs_coord_tup = {(int(tile_key[-4:]), pt[0], pt[1]): hit_given_aligned_idx[i]
                                   for i, (tile_key, pt) in enumerate(self.rcs_in_frame)}
        rcs_coord_tups = set(hit_given_rcs_coord_tup.keys())

        def flux_info_given_rcs_coord_tup(coord_tup):
            hit_type, (i, _) = hit_given_rcs_coord_tup[coord_tup]
            if hit_type == 'non_mutual':
                return 'none', 0, 0, 0, 0
            else:
                sexcat_pt = self.sexcat.points[i]
                return hit_type, sexcat_pt.r, sexcat_pt.c, sexcat_pt.flux, sexcat_pt.flux_err

        lines = set()  # set rather than list due to read pairs
        for tile in self.fastq_tiles_list:
            for read_name in tile.read_names:
                coord_tup = tuple(map(int, read_name.split(':')[-3:]))  # tile:r:c
                if coord_tup in rcs_coord_tups:
                    hit_type, rr, cc, flux, flux_err = flux_info_given_rcs_coord_tup(coord_tup)
                    lines.add('\t'.join([read_name,
                                         self.image_data.fname,
                                         hit_type,
                                         str(rr),
                                         str(cc),
                                         str(flux),
                                         str(flux_err)]))

        with open(out_fpath, 'w') as out:
            fields = ('read_name', 'image_name', 'hit_type', 'r', 'c', 'flux', 'flux_err')
            out.write('# Fields: ' + '\t'.join(fields) + '\n')
            out.write('\n'.join(sorted(lines, key=lambda s: float(s.split()[3]), reverse=True)))

    def write_alignment_stats(self, out_fpath):
        stats = [
            'Image:                 %s' % self.image_data.fname,
            'Objective:             %d' % self.image_data.objective,
            'Project Name:          %s' % self.chip_id,
            'Tile:                  %s' % ','.join(tile.key for tile in self.hitting_tiles),
            'Rotation (deg):        %s' % ','.join('%.4f' % tile.rotation_degrees for tile in self.hitting_tiles),
            'Tile width (um):       %s' % ','.join('%.4f' % tile.w for tile in self.hitting_tiles),
            'Scaling (px/fqu):      %s' % ','.join('%.7f' % tile.scale for tile in self.hitting_tiles),
            'RC Offset (px):        %s' % ','.join('(%.4f,%.4f)' % tuple(tile.offset) for tile in self.hitting_tiles),
            'Correlation:           %s' % ','.join('%.2f' % tile.best_max_corr for tile in self.hitting_tiles),
            'Non-mutual hits:       %d' % (len(self.non_mutual_hits)),
            'Bad mutual hits:       %d' % (len(self.bad_mutual_hits)),
            'Good mutual hits:      %d' % (len(self.good_mutual_hits)),
            'Exclusive hits:        %d' % (len(self.exclusive_hits)),
            'Sextractor Ellipses:   %d' % (len(self.sexcat.point_rcs)),
            'Fastq Points:          %d' % (len(self.aligned_rcs_in_frame)),
            ]
        with open(out_fpath, 'w') as out:
            out.write('\n'.join(stats))

    def write_read_names_rcs(self, out_fpath):
        im_shape = self.image_data.im.shape
        with open(out_fpath, 'w') as out:
            for tile in self.hitting_tiles:
                for read_name, pt in izip(tile.read_names, tile.aligned_rcs):
                    if 0 <= pt[0] < im_shape[0] and 0 <= pt[1] < im_shape[1]:
                        out.write('%s\t%f\t%f\n' % (read_name, pt[0], pt[1]))
