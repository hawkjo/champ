import logging
import time
from copy import deepcopy
from itertools import izip
import numpy as np
from champ import stats, clusters
from fastqtilercs import FastqTileRCs
from imagedata import ImageData
from scipy.spatial import KDTree

log = logging.getLogger(__name__)


class FastqImageAligner(object):
    """A class to find the alignment of fastq data and image data."""
    def __init__(self, microns_per_pixel):
        self.fastq_tiles = {}
        self.fastq_tiles_list = []
        self.fastq_tiles_keys = []
        self.microns_per_pixel = microns_per_pixel
        self.image_data = None
        self.fq_w = 935  # um
        self.control_corr = 0
        self.non_mutual_hits = set()
        self.mutual_hits = set()
        self.bad_mutual_hits = set()
        self.good_mutual_hits = set()
        self.exclusive_hits = set()
        self.hitting_tiles = []

    def load_reads(self, tile_data, valid_keys=None):
        for tile_key, read_names in tile_data.items():
            if valid_keys is None or tile_key in valid_keys:
                self.fastq_tiles[tile_key] = FastqTileRCs(tile_key, read_names, self.microns_per_pixel)
        self.fastq_tiles_list = [tile for tile_key, tile in sorted(self.fastq_tiles.items())]

    def all_reads_fic_from_aligned_fic(self, other_fic, all_reads):
        self.load_reads(all_reads, valid_keys=[tile.key for tile in other_fic.hitting_tiles])
        self.image_data = deepcopy(other_fic.image_data)
        self.fq_w = other_fic.fq_w
        self.set_fastq_tile_mappings()
        self.set_all_fastq_image_data()
        self.hitting_tiles = [self.fastq_tiles[tile.key] for tile in other_fic.hitting_tiles]
        self.clusters = other_fic.clusters
        for other_tile in other_fic.hitting_tiles:
            tile = self.fastq_tiles[other_tile.key]
            tile.set_aligned_rcs_given_transform(other_tile.scale,
                                                 other_tile.rotation,
                                                 other_tile.offset)

    def set_tile_alignment(self, tile_key, scale, fq_w, rotation, rc_offset):
        if self.fastq_tiles[tile_key] not in self.hitting_tiles:
            self.hitting_tiles.append(self.fastq_tiles[tile_key])
        self.fq_w = fq_w
        self.set_fastq_tile_mappings()
        self.set_all_fastq_image_data()
        tile = self.fastq_tiles[tile_key]
        tile.set_aligned_rcs_given_transform(scale, rotation, rc_offset)

    def alignment_from_alignment_file(self, path):
        self.hitting_tiles = []
        with open(path) as f:
            astats = stats.AlignmentStats().from_file(f)
        for tile_key, scaling, tile_width, rotation, rc_offset, _ in astats:
            self.set_tile_alignment(tile_key, scaling, tile_width, rotation, rc_offset)

    def set_clusters_from_file(self, fpath):
        with open(fpath) as f:
            self.clusters = clusters.Clusters(f)

    def set_image_data(self, image, um_per_pixel):
        self.image_data = ImageData(image.index, um_per_pixel, image)

    def set_all_fastq_image_data(self):
        for key, tile in self.fastq_tiles.items():
            tile.set_fastq_image_data(self.fq_im_offset,
                                      self.fq_im_scale,
                                      self.fq_im_scaled_dims,
                                      self.fq_w)

    def rotate_all_fastq_data(self, degrees):
        im_shapes = [tile.rotate_data(degrees) for tile in self.fastq_tiles_list]
        self.fq_im_scaled_dims = np.array(im_shapes).max(axis=0)
        for tile in self.fastq_tiles_list:
            tile.image_shape = self.fq_im_scaled_dims

    def set_fastq_tile_mappings(self):
        """Calculate parameters for mapping fastq tiles for ffts."""
        assert self.image_data is not None, 'No image data loaded.'
        assert self.fastq_tiles != {}, 'No fastq data loaded.'

        self.all_data = np.concatenate([tile.rcs for key, tile in self.fastq_tiles.items()])

        x_min, y_min = self.all_data.min(axis=0)
        x_max, y_max = self.all_data.max(axis=0)

        self.fq_im_offset = np.array([-x_min, -y_min])
        self.fq_im_scale = (float(self.fq_w) / (x_max-x_min)) / self.image_data.um_per_pixel
        self.fq_im_scaled_maxes = self.fq_im_scale * np.array([x_max-x_min, y_max-y_min])
        self.fq_im_scaled_dims = (self.fq_im_scaled_maxes + [1, 1]).astype(np.int)

    def find_hitting_tiles(self, possible_tile_keys, snr_thresh=1.2):
        possible_tiles = [self.fastq_tiles[key] for key in possible_tile_keys
                          if key in self.fastq_tiles]
        impossible_tiles = [tile for tile in self.fastq_tiles.values() if tile not in possible_tiles]
        impossible_tiles.sort(key=lambda tile: -len(tile.read_names))
        control_tiles = impossible_tiles[:2]
        self.image_data.set_fft(self.fq_im_scaled_dims)
        self.control_corr = 0

        for control_tile in control_tiles:
            corr, _ = control_tile.fft_align_with_im(self.image_data)
            if corr > self.control_corr:
                self.control_corr = corr
        del control_tiles
        self.hitting_tiles = []
        for tile in possible_tiles:
            max_corr, align_tr = tile.fft_align_with_im(self.image_data)
            if max_corr > snr_thresh * self.control_corr:
                tile.set_aligned_rcs(align_tr)
                tile.snr = max_corr / self.control_corr
                self.hitting_tiles.append(tile)

    def find_points_in_frame(self, consider_tiles='all'):
        self.rcs_in_frame = []
        aligned_rcs_in_frame = []
        im_shape = self.image_data.image.shape

        if consider_tiles == 'all':
            considered_tiles = self.hitting_tiles
        else:
            considered_tiles = [consider_tiles]

        for tile in considered_tiles:
            rcs = tile.rcs.astype(np.int)
            for i, pt in enumerate(tile.aligned_rcs):
                if 0 <= pt[0] < im_shape[0] and 0 <= pt[1] < im_shape[1]:
                    aligned_rcs_in_frame.append(pt)
                    self.rcs_in_frame.append((tile.key, rcs[i]))
        self.aligned_rcs_in_frame = np.array(aligned_rcs_in_frame)

    def hit_dists(self, hits):
        return [self.single_hit_dist(hit) for hit in hits]

    def single_hit_dist(self, hit):
        return np.linalg.norm(self.clusters.point_rcs[hit[0]] - self.aligned_rcs_in_frame[hit[1]])

    def remove_longest_hits(self, hits, pct_thresh):
        if not hits:
            return []
        dists = self.hit_dists(hits)
        thresh = np.percentile(dists, pct_thresh * 100)
        return [hit for hit in hits if self.single_hit_dist(hit) <= thresh]

    def find_hits(self, consider_tiles='all'):
        # --------------------------------------------------------------------------------
        # Find nearest neighbors
        # --------------------------------------------------------------------------------
        self.find_points_in_frame(consider_tiles)
        cluster_tree = KDTree(self.clusters.point_rcs)
        aligned_tree = KDTree(self.aligned_rcs_in_frame)

        # All indices are in the order (cluster_index, aligned_in_frame_idx)
        cluster_to_aligned_indexes = set()
        for i, pt in enumerate(self.clusters.point_rcs):
            dist, idx = aligned_tree.query(pt)
            cluster_to_aligned_indexes.add((i, idx))

        aligned_to_cluster_indexes_rev = set()
        for i, pt in enumerate(self.aligned_rcs_in_frame):
            dist, idx = cluster_tree.query(pt)
            aligned_to_cluster_indexes_rev.add((idx, i))

        # --------------------------------------------------------------------------------
        # Find categories of hits
        # --------------------------------------------------------------------------------
        mutual_hits = cluster_to_aligned_indexes & aligned_to_cluster_indexes_rev
        non_mutual_hits = cluster_to_aligned_indexes ^ aligned_to_cluster_indexes_rev

        clusters_in_non_mutual = set(i for i, j in non_mutual_hits)
        aligned_in_non_mutual = set(j for i, j in non_mutual_hits)
        exclusive_hits = set((i, j) for i, j in mutual_hits if i not in
                             clusters_in_non_mutual and j not in aligned_in_non_mutual)

        # --------------------------------------------------------------------------------
        # Recover good non-exclusive mutual hits. 
        # --------------------------------------------------------------------------------
        # If the distance to second neighbor is too close, that suggests a bad peak call combining
        # two peaks into one. Filter those out with a gaussian-mixture-model-determined threshold.
        if int(16.0 / self.image_data.um_per_pixel) == 60:
            # Value decided by observation of our data. May vary with equipment.

            good_hit_threshold = 5
        else:
            good_hit_threshold = np.percentile(self.hit_dists(exclusive_hits), 95)
        second_neighbor_thresh = 2 * good_hit_threshold

        exclusive_hits = set(hit for hit in exclusive_hits
                             if self.single_hit_dist(hit) <= good_hit_threshold)

        good_mutual_hits = set()
        for i, j in (mutual_hits - exclusive_hits):
            if self.hit_dists([(i, j)])[0] > good_hit_threshold:
                continue
            third_wheels = [tup for tup in non_mutual_hits if i == tup[0] or j == tup[1]]
            if min(self.hit_dists(third_wheels)) > second_neighbor_thresh:
                good_mutual_hits.add((i, j))
        bad_mutual_hits = mutual_hits - exclusive_hits - good_mutual_hits

        # --------------------------------------------------------------------------------
        # Test that the four groups form a partition of all hits and finalize
        # --------------------------------------------------------------------------------
        assert (non_mutual_hits | bad_mutual_hits | good_mutual_hits | exclusive_hits
                == cluster_to_aligned_indexes | aligned_to_cluster_indexes_rev
                and len(non_mutual_hits) + len(bad_mutual_hits)
                + len(good_mutual_hits) + len(exclusive_hits)
                == len(cluster_to_aligned_indexes | aligned_to_cluster_indexes_rev))

        self.non_mutual_hits = non_mutual_hits
        self.mutual_hits = mutual_hits
        self.bad_mutual_hits = bad_mutual_hits
        self.good_mutual_hits = good_mutual_hits
        self.exclusive_hits = exclusive_hits

        if consider_tiles != 'all':
            log.debug('Non-mutual hits: %s' % len(non_mutual_hits))
            log.debug('Mutual hits: %s' % len(mutual_hits))
            log.debug('Bad mutual hits: %s' % len(bad_mutual_hits))
            log.debug('Good mutual hits: %s' % len(good_mutual_hits))
            log.debug('Exclusive hits: %s' % len(exclusive_hits))

    def least_squares_mapping(self, pct_thresh=0.9, min_hits=50):
        """least_squares_mapping(self, hit_type='exclusive')

        "Input": set of tuples of (cluster_index, in_frame_idx) mappings.

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

        The r and s subscripts indicate rcs and cluster coords.

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

        found_good_mapping = False

        for tile in self.hitting_tiles:
            self.find_hits(consider_tiles=tile)
            # Reminder: All indices are in the order (cluster_index, in_frame_idx)
            raw_hits = get_hits(('exclusive', 'good_mutual'))
            hits = self.remove_longest_hits(raw_hits, pct_thresh)
            if len(hits) < min_hits:
                continue
            else:
                found_good_mapping = True
            A = np.zeros((2 * len(hits), 4))
            b = np.zeros((2 * len(hits),))
            for i, (cluster_index, in_frame_idx) in enumerate(hits):
                tile_key, (xir, yir) = self.rcs_in_frame[in_frame_idx]
                A[2*i, :] = [xir, -yir, 1, 0]
                A[2*i+1, :] = [yir,  xir, 0, 1]

                xis, yis = self.clusters.point_rcs[cluster_index]
                b[2*i] = xis
                b[2*i+1] = yis

            alpha, beta, x_offset, y_offset = np.linalg.lstsq(A, b)[0]
            offset = np.array([x_offset, y_offset])
            theta = np.arctan2(beta, alpha)
            lbda = alpha / np.cos(theta)
            tile.set_aligned_rcs_given_transform(lbda, theta, offset)
            tile.set_correlation(self.image_data.image)
            if hasattr(self, 'control_corr'):
                tile.set_snr_with_control_corr(self.control_corr)
        return found_good_mapping

    def rough_align(self, possible_tile_keys, rotation_est, fq_w_est=927, snr_thresh=1.2):
        self.fq_w = fq_w_est
        self.set_fastq_tile_mappings()
        self.set_all_fastq_image_data()
        self.rotate_all_fastq_data(rotation_est)
        start_time = time.time()
        self.find_hitting_tiles(possible_tile_keys, snr_thresh)
        log.debug('Rough alignment time: %.3f seconds' % (time.time() - start_time))

    def precision_align_only(self, min_hits):
        start_time = time.time()
        if not self.hitting_tiles:
            raise RuntimeError('Alignment not found')
        found_good_mapping = self.least_squares_mapping(min_hits=min_hits)
        if not found_good_mapping:
            raise ValueError("Could not precision align!")
        log.debug('Precision alignment time: %.3f seconds' % (time.time() - start_time))
        start_time = time.time()
        self.find_hits()
        log.debug('Hit finding time: %.3f seconds' % (time.time() - start_time))

    @property
    def alignment_stats(self):
        hits = {'exclusive': len(self.exclusive_hits),
                'good_mutual': len(self.good_mutual_hits),
                'bad_mutual': len(self.bad_mutual_hits),
                'non_mutual': len(self.non_mutual_hits)}
        offsets = [tuple(map(float, tile.offset)) for tile in self.hitting_tiles]
        return stats.AlignmentStats().from_data([str(tile.key) for tile in self.hitting_tiles],
                                                [float(tile.scale) for tile in self.hitting_tiles],
                                                [float(tile.width) for tile in self.hitting_tiles],
                                                [float(tile.rotation_degrees) for tile in self.hitting_tiles],
                                                offsets,
                                                hits)

    @property
    def read_names_rcs(self):
        im_shape = self.image_data.image.shape
        for tile in self.hitting_tiles:
            for read_name, pt in izip(tile.read_names, tile.aligned_rcs):
                if 0 <= pt[0] < im_shape[0] and 0 <= pt[1] < im_shape[1]:
                    yield '%s\t%f\t%f\n' % (read_name, pt[0], pt[1])
