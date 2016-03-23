from chimp.model.sextractor import Sextraction
from chimp.model.grid import GridImages
from chimp.model import tiles
from chimp import fastq
from collections import defaultdict
from functools import partial
import logging
from nd2reader import Nd2
import numpy as np
import os
import padding
import random
import scipy.optimize
from scipy.spatial import KDTree
from sklearn.mixture import GMM
from skimage import io
import matplotlib.pyplot as plt


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
log.addHandler(logging.StreamHandler())


def load_sexcat(directory, index):
    with open(os.path.join(directory, '%s.cat' % index)) as f:
        return Sextraction(f)


def max_2d_idx(a):
    return np.unravel_index(a.argmax(), a.shape)


def correlate_image_and_tile(image, tile_fft_conjugate):
    dimension = padding.calculate_pad_size(tile_fft_conjugate.shape[0],
                                           tile_fft_conjugate.shape[1],
                                           image.shape[0],
                                           image.shape[1])
    padded_microscope = padding.pad_image(image, dimension)
    image_fft = np.fft.fft2(padded_microscope)
    cross_correlation = abs(np.fft.ifft2(tile_fft_conjugate * image_fft))
    max_index = max_2d_idx(cross_correlation)
    alignment_transform = np.array(max_index) - image_fft.shape
    return cross_correlation.max(), max_index, alignment_transform


def get_expected_tile_map(min_tile, max_tile, min_column, max_column):
    """
    Creates a dictionary that relates each column of microscope images to its expected tile, +/- 1.

    """
    # Mi-Seq chips have 19 tiles on each side
    MISEQ_TILE_COUNT = 19
    tile_map = defaultdict(list)
    normalization_factor = float(max_tile - min_tile + 1) / float(max_column - min_column)
    for column in range(min_column, max_column + 1):
        expected_tile = min(MISEQ_TILE_COUNT, max(1, int(round(column * normalization_factor, 0)))) + min_tile - 1
        tile_map[column].append(expected_tile)
        if expected_tile > min_tile:
            tile_map[column].append(expected_tile - 1)
        if expected_tile < max_tile:
            tile_map[column].append(expected_tile + 1)
    return tile_map


def main(base_image_name, snr, alignment_channel=None, alignment_offset=None, cache_size=20):
    LEFT_TILES = tuple(range(1, 11))
    RIGHT_TILES = tuple(reversed(range(11, 20)))
    nd2 = Nd2('%s.nd2' % base_image_name)
    loader = partial(load_sexcat, base_image_name)
    log.debug("Loading mapped reads.")
    mapped_reads = fastq.load_mapped_reads('phix')
    tm = tiles.load_tile_manager(nd2.pixel_microns, nd2.height, nd2.width, mapped_reads, cache_size)
    grid = GridImages(nd2, loader, alignment_channel, alignment_offset)
    log.debug("Loaded grid and tile manager.")
    left_tile, left_column = find_end_tile(grid.left_iter(), tm, snr, LEFT_TILES, RIGHT_TILES, "left")
    right_tile, right_column = find_end_tile(grid.right_iter(), tm, snr, RIGHT_TILES, LEFT_TILES, "right")
    log.debug("Data was acquired between tiles %s and %s" % (left_tile, right_tile))
    log.debug("Using images from column %s and %s" % (left_column, right_column))

    # get a dictionary of tiles that each image will probably be found in
    tile_map = get_expected_tile_map(left_tile, right_tile, left_column, right_column)

    # results is a hacky temporary data structure for storing alignment results. Jim WILL make it better
    results = defaultdict(list)
    for microscope_data in grid.bounded_iter(left_column, right_column):
        for tile_number, cross_correlation, max_index, alignment_transform in get_rough_alignment(microscope_data, tm, tile_map, snr):
            results[(microscope_data.row, microscope_data.column)].append((tile_number, cross_correlation, max_index, alignment_transform))
            tile = tm.get(tile_number)
            hits = find_hits(microscope_data, tile, alignment_transform)


def find_hits(microscope_data, tile, rough_alignment_transform):
    tile_points = []
    point_indexes = []
    for i, original_point in enumerate(tile.normalized_rcs):
        point = original_point + rough_alignment_transform
        if 0 <= point[0] < microscope_data.shape[0] and 0 <= point[1] < microscope_data.shape[1]:
            tile_points.append(point)
            point_indexes.append(i)

    aligned_rcs = np.array(tile_points).astype(np.int)
    sexcat_tree = KDTree(microscope_data.rcs)
    aligned_tree = KDTree(aligned_rcs)
    sexcat_to_aligned_idxs = set()
    for i, pt in enumerate(microscope_data.rcs):
        dist, idx = aligned_tree.query(pt)
        sexcat_to_aligned_idxs.add((i, idx))

    aligned_to_sexcat_idxs_rev = set()
    for i, pt in enumerate(aligned_rcs):
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
    good_hit_thresh = 5
    second_neighbor_thresh = 2 * good_hit_thresh
    exclusive_hits = set(hit for hit in exclusive_hits
                         if single_hit_dist(hit, microscope_data.rcs, aligned_rcs) <= good_hit_thresh)

    good_mutual_hits = set()
    for i, j in (mutual_hits - exclusive_hits):
        if single_hit_dist((i, j), microscope_data.rcs, aligned_rcs) > good_hit_thresh:
            continue
        third_wheels = [tup for tup in non_mutual_hits if i == tup[0] or j == tup[1]]
        third_wheel_distances = [single_hit_dist(hit, microscope_data.rcs, aligned_rcs) for hit in third_wheels]
        if min(third_wheel_distances) > second_neighbor_thresh:
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

    log.debug('Non-mutual hits: %s' % len(non_mutual_hits))
    log.debug('Mutual hits: %s' % len(mutual_hits))
    log.debug('Bad mutual hits: %s' % len(bad_mutual_hits))
    log.debug('Good mutual hits: %s' % len(good_mutual_hits))
    log.debug('Exclusive hits: %s' % len(exclusive_hits))
    return exclusive_hits


def single_hit_dist(hit, microscope_rcs, aligned_rcs):
    return np.linalg.norm(microscope_rcs[hit[0]] - aligned_rcs[hit[1]])


def get_rough_alignment(microscope_data, tile_manager, tile_map, snr):
    possible_tiles = set(tile_map[microscope_data.column])
    impossible_tiles = set(range(1, 20)) - possible_tiles
    control_correlation = calculate_control_correlation(microscope_data.image, tile_manager, impossible_tiles)
    noise_threshold = control_correlation * snr
    for tile_number in possible_tiles:
        cross_correlation, max_index, alignment_transform = correlate_image_and_tile(microscope_data.image,
                                                                                     tile_manager.fft_conjugate(tile_number))
        if cross_correlation > noise_threshold:
            log.debug("Field of view %sx%s is in tile %s" % (microscope_data.row, microscope_data.column, tile_number))
            yield tile_number, cross_correlation, max_index, alignment_transform


def calculate_control_correlation(image, tile_manager, impossible_tiles):
    control_correlation = 0.0
    for impossible_tile in random.sample(impossible_tiles, 3):
        tile_fft_conjugate = tile_manager.fft_conjugate(impossible_tile)
        cross_correlation, max_index, alignment_transform = correlate_image_and_tile(image, tile_fft_conjugate)
        control_correlation = max(control_correlation, cross_correlation)
    return control_correlation


def find_end_tile(grid_images, tile_manager, snr, possible_tiles, impossible_tiles, side_name):
    log.debug("Finding end tiles on the %s" % side_name)
    for microscope_data in grid_images:
        # first get the correlation to random tiles, so we can distinguish signal from noise
        control_correlation = calculate_control_correlation(microscope_data.image, tile_manager, impossible_tiles)
        # now find which tile the image aligns to, if any
        for tile_number in possible_tiles:
            tile_fft_conjugate = tile_manager.fft_conjugate(tile_number)
            cross_correlation, max_index, alignment_transform = correlate_image_and_tile(microscope_data.image,
                                                                                         tile_fft_conjugate)
            if cross_correlation > (control_correlation * snr):
                return tile_number, microscope_data.column


if __name__ == '__main__':
    results = main('/var/experiments/151118/15-11-18_SA15243_Cascade-TA_1nM-007', 1.2, alignment_offset=1, cache_size=8)
