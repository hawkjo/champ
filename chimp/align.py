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
from scipy.spatial import KDTree

log = logging.getLogger()
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


def main(base_image_name, snr, alignment_channel=None, alignment_offset=None, cache_size=20,
         precision_hit_threshold=0.9, lane=1, side=2):
    log.info("Starting alignment process for %s" % base_image_name)
    LEFT_TILES = tuple(range(1, 11))
    RIGHT_TILES = tuple(reversed(range(11, 20)))
    nd2 = Nd2('%s.nd2' % base_image_name)
    loader = partial(load_sexcat, base_image_name)
    phix_reads = fastq.load_mapped_reads('phix')
    tm = tiles.load_tile_manager(nd2.pixel_microns, nd2.height, nd2.width, phix_reads, cache_size)
    grid = GridImages(nd2, loader, alignment_channel, alignment_offset)
    # left_tile, left_column = find_end_tile(grid.left_iter(), tm, snr, LEFT_TILES, RIGHT_TILES, "left")
    # right_tile, right_column = find_end_tile(grid.right_iter(), tm, snr, RIGHT_TILES, LEFT_TILES, "right")
    left_tile, left_column, right_tile, right_column = 4, 0, 14, 59
    log.debug("Data was acquired between tiles %s and %s" % (left_tile, right_tile))
    log.debug("Using images from column %s and %s" % (left_column, right_column))

    # get a dictionary of tiles that each image will probably be found in
    tile_map = get_expected_tile_map(left_tile, right_tile, left_column, right_column)
    images = grid.bounded_iter(left_column, right_column)

    results = align(images, tm, tile_map, snr, precision_hit_threshold, lane, side)

    for index, rs in results.items():
        for result in rs:
            fastq_read_rcs = {}
            for read in phix_reads.get((result.lane, result.side, result.tile_number)):
                fastq_read_rcs[(int(read.column), int(read.row))] = read.name

            tile = tm.get(result.tile_number)
            corrected_rcs = tiles.flip_coordinates(result.original_rcs)
            corrected_rcs = tiles.normalize_rcs(tile.scale, tile.offset, corrected_rcs)
            corrected_rcs -= (tile.scale * tile.offset)
            # put next line after precision align?
            corrected_rcs = tiles.precision_align_rcs(corrected_rcs,
                                                      result.offset,
                                                      result.theta,
                                                      result.lbda)
            # for i, ((o_r, o_c), (c_r, c_c)) in enumerate(zip(result.original_rcs, corrected_rcs)):
            #     name = fastq_read_rcs.get((o_r, o_c))
            #     if name == "M02288:175:000000000-AHFHH:1:2104:3490:17288":
            #         print(i, name, o_r, o_c, c_r / 452.527219, c_c / 188.160377)


class Result(object):
    def __init__(self, index, row, column, lane, side, tile_number):
        self.index = index
        self.row = row
        self.column = column
        self.lane = lane
        self.side = side
        self.tile_number = tile_number
        self.offset = None
        self.theta = None
        self.at = None
        self.lbda = None
        self.original_rcs = None

    def set_precision_alignment(self, offset, theta, lbda):
        self.offset = offset
        self.theta = theta
        self.lbda = lbda

    def set_rcs(self, original_rcs):
        self.original_rcs = original_rcs


def align(images, tm, tile_map, snr, precision_hit_threshold, lane, side):
    results = defaultdict(list)
    for n, microscope_data in enumerate(images):
        log.debug("aligning %sx%s" % (microscope_data.row, microscope_data.column))
        for tile_number, cross_correlation, max_index, alignment_transform in get_rough_alignment(microscope_data, tm, tile_map, snr):
            # the rough alignment worked for this tile, so now calculate the precision alignment
            tile = tm.get(tile_number, lane, side)
            rough_aligned_rcs, original_rcs = find_aligned_rcs(microscope_data, tile, alignment_transform)
            exclusive_hits = find_hits(microscope_data, rough_aligned_rcs, precision_hit_threshold)
            offset, theta, lbda = least_squares_mapping(exclusive_hits, microscope_data.sexcat_rcs, rough_aligned_rcs)

            # save the results
            result = Result(microscope_data.index, microscope_data.row, microscope_data.column, lane, side, tile_number)
            result.set_rcs(original_rcs)
            result.set_precision_alignment(offset, theta, lbda)
            results[microscope_data.index].append(result)
        # skip for testing
        if n == 1:
            return results
    log.debug("Done aligning!")
    return results


def remove_longest_hits(hits, microscope_data, aligned_rcs, pct_thresh):
    dists = [single_hit_dist(hit, microscope_data.sexcat_rcs, aligned_rcs) for hit in hits]
    proximity_threshold = np.percentile(dists, pct_thresh * 100)
    nice_hits = [hit for hit in hits if single_hit_dist(hit, microscope_data.sexcat_rcs, aligned_rcs) <= proximity_threshold]
    return nice_hits


def least_squares_mapping(hits, sexcat_rcs, inframe_tile_rcs, min_hits=50):
    """
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
    # Reminder: All indices are in the order (sexcat_idx, in_frame_idx)
    assert len(hits) > min_hits, 'Too few hits for least squares mapping: {0}'.format(len(hits))
    A = np.zeros((2 * len(hits), 4))
    b = np.zeros((2 * len(hits),))
    for i, (sexcat_idx, in_frame_idx) in enumerate(hits):
        xir, yir = inframe_tile_rcs[in_frame_idx]
        A[2*i, :] = [xir, -yir, 1, 0]
        A[2*i+1, :] = [yir,  xir, 0, 1]

        xis, yis = sexcat_rcs[sexcat_idx]
        b[2*i] = xis
        b[2*i+1] = yis

    alpha, beta, x_offset, y_offset = np.linalg.lstsq(A, b)[0]
    offset = np.array([x_offset, y_offset])
    theta = np.arctan2(beta, alpha)
    lbda = alpha / np.cos(theta)
    return offset, theta, lbda


def find_aligned_rcs(microscope_data, tile, rough_alignment_transform):
    tile_points = []
    original_rcs = []
    for normalized_point, raw_point in zip(tile.normalized_rcs, tile.raw_rcs):
        point = normalized_point + rough_alignment_transform
        if 0 <= point[0] < microscope_data.shape[0] and 0 <= point[1] < microscope_data.shape[1]:
            tile_points.append(point)
            original_rcs.append(raw_point)
    return np.array(tile_points).astype(np.int64), np.array(original_rcs).astype(np.int64)


def find_hits(microscope_data, aligned_rcs, precision_hit_threshold):
    sexcat_tree = KDTree(microscope_data.sexcat_rcs)
    aligned_tree = KDTree(aligned_rcs)
    sexcat_to_aligned_idxs = set()
    for i, pt in enumerate(microscope_data.sexcat_rcs):
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

    good_hit_thresh = 5
    exclusive_hits = set(hit for hit in exclusive_hits
                         if single_hit_dist(hit, microscope_data.sexcat_rcs, aligned_rcs) <= good_hit_thresh)
    return remove_longest_hits(exclusive_hits, microscope_data, aligned_rcs, precision_hit_threshold)


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
    main('/var/experiments/151118/15-11-18_SA15243_Cascade-TA_1nM-007', 1.2, alignment_offset=1, cache_size=7)
