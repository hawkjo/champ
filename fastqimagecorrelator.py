import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import ndimage
import local_config
import sextraction
import scipy.optimize
from scipy.spatial import KDTree
from sklearn.mixture import GMM
from Bio import SeqIO
from imagedata import ImageData
from fastqtilercs import FastqTileRCs
from misc import pad_to_size, max_2d_idx


class FastqImageCorrelator(object):
    """A class to find the alignment of fastq data and image data."""
    def __init__(self, project_name):
        self.project_name = project_name
        self.fastq_tiles = {}
        self.fastq_tiles_list = []
        self.fastq_tiles_keys = []
        self.image_data = ImageData()
        self.w_fq_tile_min = 895  # um
        self.w_fq_tile_max = 937  # um
        self.w_fq_tile = 937  # um

    def load_phiX(self):
        for key, tile in local_config.phiX_rcs_given_project_name(self.project_name).items():
            self.fastq_tiles[key] = FastqTileRCs(key, tile)
        self.fastq_tiles_keys = [key for key, tile in sorted(self.fastq_tiles.items())]
        self.fastq_tiles_list = [tile for key, tile in sorted(self.fastq_tiles.items())]

    def set_image_data(self, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], ImageData):
            self.image_data = args[0]
        else:
            self.image_data = ImageData(*args, **kwargs)

    def set_fastq_tile_mappings(self):
        """Calculate parameters for mapping fastq tiles for ffts."""
        assert self.image_data is not None, 'No image data loaded.'
        assert self.fastq_tiles != {}, 'No fastq data loaded.'

        self.all_data = np.concatenate([tile.rcs for key, tile in self.fastq_tiles.items()])
    
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
        #aligned_rcs = self.fastq_tiles[fq_key].aligned_rcs()
        #ax.plot(aligned_rcs[:, 0], aligned_rcs[:, 1], 'r.', alpha=0.5)
        plt.axis(v)

    def set_sexcat(self, sexcat):
        assert isinstance(sexcat, sextraction.Sextraction)
        self.sexcat = sexcat

    def set_sexcat_from_file(self, fpath):
        self.sexcat = sextraction.Sextraction(fpath)

    def find_hitting_tiles(self, possible_tile_keys, snr_thresh=2):
        possible_tiles = [self.fastq_tiles[key] for key in possible_tile_keys]
        for tile in self.fastq_tiles_list:
            if tile not in possible_tiles:
                control_tile = tile
                break

        self.hitting_tiles = []
        self.image_data.set_single_fft((0, 0), padding=self.fq_im_scaled_dims)
        _, _, control_corr, _ = self.fft_align_tile(control_tile)
        for tile in possible_tiles:
            tile_key, _, max_corr, align_tr = self.fft_align_tile(tile)
            if max_corr > snr_thresh * control_corr:
                tile.set_aligned_rcs()
                self.hitting_tiles.append(tile)

    def find_points_in_frame(self):
        self.aligned_rcs_in_frame = []
        self.rcs_in_frame = []
        im_shape = self.image_data.im.shape
        for tile in self.hitting_tiles:
            rcs = tile.rcs.astype(np.int)
            for i, pt in enumerate(tile.aligned_rcs):
                if 0 <= pt[0] < im_shape[0] and 0 <= pt[1] < im_shape[1]:
                    self.aligned_rcs_in_frame.append(pt)
                    self.rcs_in_frame.append((tile.key, rcs[i]))
        self.aligned_rcs_in_frame = np.array(self.aligned_rcs_in_frame)

    def hit_dists(self, hits):
        dists = []
        for i, j in hits:
            dists.append(np.linalg.norm(self.sexcat.point_rcs[i] - self.aligned_rcs_in_frame[j]))
        return dists

    def find_hits(self, good_posterior_thresh_pct=0.99, good_mutual_hit_thresh=None):
        #--------------------------------------------------------------------------------
        # Find nearest neighbors
        #--------------------------------------------------------------------------------
        self.find_points_in_frame()
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

        #--------------------------------------------------------------------------------
        # Find categories of hits
        #--------------------------------------------------------------------------------
        mutual_hits = sexcat_to_aligned_idxs & aligned_to_sexcat_idxs_rev
        non_mutual_hits = sexcat_to_aligned_idxs ^ aligned_to_sexcat_idxs_rev

        sexcat_in_non_mutual = set(i for i, j in non_mutual_hits)
        aligned_in_non_mutual = set(j for i, j in non_mutual_hits)
        exclusive_hits = set((i, j) for i, j in mutual_hits if i not in
                             sexcat_in_non_mutual and j not in aligned_in_non_mutual)

        #--------------------------------------------------------------------------------
        # Recover good non-exclusive mutual hits. 
        #--------------------------------------------------------------------------------
        # If the distance to second neighbor is too close, that suggests a bad peak call combining
        # two peaks into one. Filter those out with a gaussian-mixture-model-determined threshold.
        if good_mutual_hit_thresh is not None:
            assert good_mutual_hit_thresh > 0
            self.good_mutual_hit_thresh = good_mutual_hit_thresh
        else:
            self.gmm_thresh(self.hit_dists(non_mutual_hits))

        good_mutual_hits = set()
        for i, j in (mutual_hits - exclusive_hits):
            third_wheels = [tup for tup in non_mutual_hits if i == tup[0] or j == tup[1]]
            if min(self.hit_dists(third_wheels)) > self.good_mutual_hit_thresh:
                good_mutual_hits.add((i, j))
        bad_mutual_hits = mutual_hits - exclusive_hits - good_mutual_hits

        #--------------------------------------------------------------------------------
        # Test that the four groups form a partition of all hits and finalize
        #--------------------------------------------------------------------------------
        assert (non_mutual_hits | bad_mutual_hits | good_mutual_hits | exclusive_hits
                == sexcat_to_aligned_idxs | aligned_to_sexcat_idxs_rev
                and  len(non_mutual_hits) + len(bad_mutual_hits)
                + len(good_mutual_hits) + len(exclusive_hits)
                == len(sexcat_to_aligned_idxs | aligned_to_sexcat_idxs_rev))
                                
        self.non_mutual_hits = non_mutual_hits
        self.mutual_hits = mutual_hits
        self.bad_mutual_hits = bad_mutual_hits
        self.good_mutual_hits = good_mutual_hits
        self.exclusive_hits = exclusive_hits

        print 'Non-mutual hits:', len(non_mutual_hits)
        print 'Mutual hits:', len(mutual_hits)
        print 'Bad mutual hits:', len(bad_mutual_hits)
        print 'Good mutual hits:', len(good_mutual_hits)
        print 'Exclusive hits:', len(exclusive_hits)

    def gmm_thresh(self, dists):
        self.gmm = GMM(2)
        self.gmm.fit(dists)
        
        lower_idx = self.gmm.means_.argmin()
        higher_idx = 1 - lower_idx

        lower_mean = self.gmm.means_[lower_idx]
        good_posterior_thresh_pct = 0.99
        f = lambda x: self.gmm.predict_proba([x])[0][higher_idx] - good_posterior_thresh_pct
        self.good_mutual_hit_thresh = scipy.optimize.brentq(f, lower_mean, max(dists))

    def plot_hit_hists(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 8))
        ax.hist(self.hit_dists(self.non_mutual_hits), 40, label='Non-mutual hits', normed=True, histtype='step')
        ax.hist(self.hit_dists(self.bad_mutual_hits), 40, label='Bad mutual hits', normed=True, histtype='step')
        ax.hist(self.hit_dists(self.good_mutual_hits), 40, label='Good mutual hits', normed=True, histtype='step')
        ax.hist(self.hit_dists(self.exclusive_hits), 40, label='Exclusive hits', normed=True, histtype='step')
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
        axs[0].plot([self.good_mutual_hit_thresh, self.good_mutual_hit_thresh], ylim,
                'g--', label='Threshold')
        axs[0].set_title('%s GMM PDF of Non-mutual hits' % self.image_data.bname)
        axs[0].legend()
        axs[0].set_ylim(ylim)

        axs[1].hist(non_mut_dists, 40, histtype='step', normed=True, label='Data')
        axs[1].plot(xs, posteriors, label='Posterior')
        axs[1].plot([self.good_mutual_hit_thresh, self.good_mutual_hit_thresh], [0, 1],
                'g--', label='Threshold')
        axs[1].set_title('%s GMM Posterior Probabilities' % self.image_data.bname)
        axs[1].legend()
        return axs

    def extract_intensity_and_sequence_from_fastq(self, fastq_fpath, out_fpath):
        hit_given_aligned_idx = {}
        for hit_type in ['non_mutual', 'bad_mutual', 'good_mutual', 'exclusive']:
            for i, j in getattr(self, hit_type + '_hits'):
                hit_given_aligned_idx[j] = (hit_type, (i, j))

        hit_given_rcs_coord_tup = {(int(tile_key[-4:]), pt[0], pt[1]): hit_given_aligned_idx[i]
                                   for i, (tile_key, pt) in enumerate(self.rcs_in_frame)}
        rcs_coord_tups = set(hit_given_rcs_coord_tup.keys())

        def flux_info_given_rcs_coord_tup(coord_tup):
            hit_type, (i, _) = hit_given_rcs_coord_tup[coord_tup]
            if hit_type == 'non_mutual':
                return 'none', 0, 0
            else:
                sexcat_pt = self.sexcat.points[i]
                return hit_type, sexcat_pt.flux, sexcat_pt.flux_err

        lines = set()  # set rather than list due to read pairs
        for record in SeqIO.parse(open(fastq_fpath), 'fastq'):
            coord_tup = tuple(map(int, str(record.id).split(':')[-3:]))  # tile:r:c
            if coord_tup in rcs_coord_tups:
                hit_type, flux, flux_err = flux_info_given_rcs_coord_tup(coord_tup)
                lines.add('\t'.join([record.id,
                                     self.image_data.fname,
                                     hit_type,
                                     str(flux),
                                     str(flux_err)]))

        with open(out_fpath, 'w') as out:
            fields = ['read_name', 'image_name', 'hit_type', 'flux', 'flux_err']
            out.write('# Fields: ' + '\t'.join(fields) + '\n')
            out.write('\n'.join(lines))

    def extract_intensity_and_sequence_from_name_list(self, name_list_fpath, out_fpath):
        hit_given_aligned_idx = {}
        for hit_type in ['non_mutual', 'bad_mutual', 'good_mutual', 'exclusive']:
            for i, j in getattr(self, hit_type + '_hits'):
                hit_given_aligned_idx[j] = (hit_type, (i, j))

        hit_given_rcs_coord_tup = {(int(tile_key[-4:]), pt[0], pt[1]): hit_given_aligned_idx[i]
                                   for i, (tile_key, pt) in enumerate(self.rcs_in_frame)}
        rcs_coord_tups = set(hit_given_rcs_coord_tup.keys())

        def flux_info_given_rcs_coord_tup(coord_tup):
            hit_type, (i, _) = hit_given_rcs_coord_tup[coord_tup]
            if hit_type == 'non_mutual':
                return 'none', 0, 0
            else:
                sexcat_pt = self.sexcat.points[i]
                return hit_type, sexcat_pt.flux, sexcat_pt.flux_err

        lines = set()  # set rather than list due to read pairs
        for line in open(name_list_fpath):
            coord_tup = tuple(map(int, line.strip().split(':')[-3:]))  # tile:r:c
            if coord_tup in rcs_coord_tups:
                hit_type, flux, flux_err = flux_info_given_rcs_coord_tup(coord_tup)
                lines.add('\t'.join([line.strip(),
                                     self.image_data.fname,
                                     hit_type,
                                     str(flux),
                                     str(flux_err)]))

        with open(out_fpath, 'w') as out:
            fields = ['read_name', 'image_name', 'hit_type', 'flux', 'flux_err']
            out.write('# Fields: ' + '\t'.join(fields) + '\n')
            out.write('\n'.join(lines))

    def plot_hits(self, hits, color, ax):
        for i, j in hits:
            ax.plot([self.sexcat.point_rcs[i, 1], self.aligned_rcs_in_frame[j, 1]],
                    [self.sexcat.point_rcs[i, 0], self.aligned_rcs_in_frame[j, 0]],
                    color=color)

    def plot_all_hits(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(15, 15))
        ax.imshow(self.image_data.im, cmap=plt.get_cmap('Blues'))
        ax.plot(self.aligned_rcs_in_frame[:, 1], self.aligned_rcs_in_frame[:, 0],
                'ko', alpha=0.3, linewidth=0, markersize=3)
        self.sexcat.plot_ellipses(ax=ax, alpha=0.6, color='darkgoldenrod')
        self.plot_hits(self.non_mutual_hits, 'grey', ax)
        self.plot_hits(self.bad_mutual_hits, 'b', ax)
        self.plot_hits(self.good_mutual_hits, 'magenta', ax)
        self.plot_hits(self.exclusive_hits, 'r', ax)
        ax.set_title('All Hits: %s vs. %s'
                % (self.image_data.bname, ','.join(tile.key for tile in self.hitting_tiles)))
        ax.set_xlim([0, self.image_data.im.shape[1]])
        ax.set_ylim([self.image_data.im.shape[0], 0])
        
        grey_line = Line2D([], [], color='grey',
                label='Non-mutual hits: %d' % (len(self.non_mutual_hits)))
        blue_line = Line2D([], [], color='blue',
                label='Bad mutual hits: %d' % (len(self.bad_mutual_hits)))
        magenta_line = Line2D([], [], color='magenta',
                label='Good mutual hits: %d' % (len(self.good_mutual_hits)))
        red_line = Line2D([], [], color='red',
                label='Exclusive hits: %d' % (len(self.exclusive_hits)))
        sexcat_line = Line2D([], [], color='darkgoldenrod', alpha=0.6, marker='o', markersize=10,
                label='Sextractor Ellipses: %d' % (len(self.sexcat.point_rcs)))
        fastq_line = Line2D([], [], color='k', alpha=0.3, marker='o', markersize=10,
                label='Fastq Points: %d' % (len(self.aligned_rcs_in_frame)))
        handles = [grey_line, blue_line, magenta_line, red_line, sexcat_line, fastq_line]
        legend = ax.legend(handles=handles)
        legend.get_frame().set_color('white')
        return ax
