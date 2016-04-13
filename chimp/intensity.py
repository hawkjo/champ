import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import nd2reader
import os
import misc
from collections import defaultdict, Counter
import time


def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    from scipy.spatial import Delaunay
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0


def alignment_obviously_bad(stats_fpath, rcs_fpath):
    """
    Obviously bad if has 3+ tiles, non-neighboring tiles, or two neighboring tiles with overlapping
    aligned point clouds.
    """
    align_stats = misc.AlignmentStats(stats_fpath)
    if align_stats.numtiles == 1:
        return False
    if align_stats.numtiles > 2:
        return True
    tile_nums = [int(tile_name[-4:]) for tile_name in align_stats.tiles]
    if abs(tile_nums[0] - tile_nums[1]) != 1:
        return True

    read_names, points = misc.read_names_and_points_given_rcs_fpath(rcs_fpath)
    grouped_points = defaultdict(list)
    for read_name, point in zip(read_names, points):
        tile_num = read_name.split(':')[4][-4:]
        grouped_points[tile_num].append(point)
    assert len(grouped_points.keys()) == 2, grouped_points.keys()
    return in_hull(*grouped_points.values())


def sum_of_gaussians(points, im_shape, params):
    """
    Calculate sum of gaussians.

    `points` is the matrix of center points for the gaussians.
    `im_shape` is the shape of the grid on which to calculate function values.
    `params` is a 1 + 2*num_points length vector containing background intensity and sigma and
        amplitude values for the gaussian at each point.
    """
    assert len(params) == 1 + 2 * points.shape[0], (len(params), points.shape[0])
    bg = params[0]
    out = bg * np.ones(im_shape)
    i = 1
    for point in points:
        r0, c0 = point
        sigma = params[i]
        A = params[i+1]
        i += 2

        r, c = np.arange(im_shape[0]), np.arange(im_shape[1])
        gr = np.exp(-(r-r0)**2/(2*sigma**2))
        gc = np.exp(-(c-c0)**2/(2*sigma**2))
        out += A * np.outer(gr, gc)
    return out


def fit_saturation_tolerant_sum_of_gaussians(im, points, sat_val=None, alpha=0.001):
    if sat_val is None:
        sat_val = im.max()
    im_eq_sat_val = im == sat_val

    def saturation_tolerant_gaussian_sq_err(params, im, points, alpha):
        """
        Calculate square error of sum of gaussians vs. pixel values, except where the sum of
        gaussians is larger than saturated pixels.

        `params` is that needed by sum_of_gaussians.
        """
        diff = im - sum_of_gaussians(points, im.shape, params)
        diff[np.logical_and(im_eq_sat_val, diff > 0)] = 0
        return (diff**2).sum() + alpha * sum(abs(params))

    x0 = [1]*(1 + 2 * points.shape[0])
    # Force sigma to be positive, bg and amplitude to be non-negative
    bounds = [(0, None) if i % 2 == 0 else (0.00001, None) for i in range(len(x0))]
    return scipy.optimize.minimize(saturation_tolerant_gaussian_sq_err,
                                   x0,
                                   args=(im, points, alpha),
                                   bounds=bounds)


def areas_given_params(params):
    """
    Returns areas of gaussians in sum_of_gaussians.
    """
    areas = []
    for i in range(1, len(params), 2):
        sigma = params[i]
        A = params[i+1]
        areas.append(2*np.pi*A*sigma**2)
    return np.array(areas)


def bin_read_names_and_points(read_names, points, im_shape, usable_side=22, boundary=3):
    """
    Bins read names and points into overlapping bins.

    `usable_side` is the side length of the usable area of each bin.
    `boundary` is the number of pixels added to each side of the usable area as a discarded
        boundary.
    """
    # Make bins
    rc_bins = []
    for idx in [0, 1]:
        edges = range(boundary, im_shape[idx], usable_side)
        if im_shape[idx] - edges[-1] < usable_side:
            # Prefer to have bins larger than `usable_side` rather than smaller.
            edges.pop()
        edges.append(im_shape[idx])
        rc_bins.append(
            [(max(e1-boundary, 0), min(e2+boundary, im_shape[idx]))
             for e1, e2 in zip(edges, edges[1:])]
        )

    # Bin
    binned_read_names, binned_adjusted_points = defaultdict(list), defaultdict(list)
    for read_name, point in zip(read_names, points):
        for rbin in rc_bins[0]:
            for cbin in rc_bins[1]:
                if rbin[0] <= point[0] < rbin[1] and cbin[0] <= point[1] < cbin[1]:
                    binned_read_names[(rbin, cbin)].append(read_name)
                    binned_adjusted_points[(rbin, cbin)].append(point)

    # Adjust points to map (rbin[0], cbin[0]) -> (0, 0)
    for key in binned_adjusted_points.keys():
        rbin, cbin = key
        bin_offset = np.array((rbin[0], cbin[0]))
        binned_adjusted_points[key] = np.array(binned_adjusted_points[key]) - bin_offset

    return binned_read_names, binned_adjusted_points


def infer_nd2_im_point_areas(nd2_fpath, im_idx, results_dir, fig_dir):
    """
    Infers all point areas in given nd2_fpath file at im_idx.
    """
    start_time = time.time()
    nd2 = nd2reader.Nd2(nd2_fpath)

    # Check for good alignment
    print 'Check alignment'
    stats_fpath = os.path.join(results_dir, '%d_stats.txt' % im_idx)
    rcs_fpath = os.path.join(results_dir, '%d_all_read_rcs.txt' % im_idx)
    if (not os.path.isfile(stats_fpath)) or alignment_obviously_bad(stats_fpath, rcs_fpath):
        return

    # Bin points
    print 'Binning points'
    im = nd2[im_idx].data
    read_names, points = misc.read_names_and_points_given_rcs_fpath(rcs_fpath)
    usable_side_len = 22
    boundary = 3
    binned_read_names, binned_adjusted_points = bin_read_names_and_points(read_names,
                                                                          points,
                                                                          im.shape,
                                                                          usable_side_len,
                                                                          boundary)

    # Infer intensities
    print 'Infer intensities'
    sigmas, amplitudes, areas = {}, {}, {}
    stats = Counter()
    for (rbin, cbin), bin_read_names in sorted(binned_read_names.items()):
        print rbin, cbin, len(bin_read_names),
        bin_points = binned_adjusted_points[(rbin, cbin)]
        bin_im = im[rbin[0]:rbin[1], cbin[0]:cbin[1]]
        res = fit_saturation_tolerant_sum_of_gaussians(bin_im, bin_points)
        if not res.success:
            stats['Failed'] += 1
            print 'Failed'
        else:
            stats['Success'] += 1
            print 'Success'
        bin_areas = areas_given_params(res.x)
        assert len(bin_read_names) == len(bin_points) == len(bin_areas) > 0, (len(bin_read_names), len(bin_points), len(bin_areas))
        rmax = rbin[1] - rbin[0] - boundary
        cmax = cbin[1] - cbin[0] - boundary
        for i, (read_name, point, area) in enumerate(zip(bin_read_names, bin_points, bin_areas)):
            if boundary <= point[0] < rmax and boundary <= point[1] < cmax:  # Discard boundary
                sigmas[read_name] = res.x[i+1]
                amplitudes[read_name] = res.x[i+2]
                areas[read_name] = area
    for stat, count in stats.items():
        print 'Total %s: %d' % (stat, count)
    print 'Inferrence time:', time.time() - start_time

    # Write areas
    print 'Write areas'
    areas_fpath = os.path.join(results_dir, '%d_gaussian_areas.txt' % im_idx)
    with open(areas_fpath, 'w') as out:
        out.write('# read_name\tsigma\tamplitude\tarea\n')
        out.write('\n'.join(['%s\t%g\t%g\t%g' % (rn, sigmas[rn], amplitudes[rn], areas[rn])
                             for rn in areas.keys()]))

    # Compare against sextractor
    print 'Compare against sextractor'
    sex_intensity_fpath = os.path.join(results_dir, '%d_intensities.txt' % im_idx)
    fig_fpath = os.path.join(fig_dir, '%d_gaussian_areas_vs_sextractor_fluxes.png' % im_idx)
    area_given_read_name = {read_name: area for read_name, area in zip(read_names, areas)}
    sex_intensities_given_read_name = {}
    for line in open(sex_intensity_fpath):
        read_name, im_name, hit_type, r, c, flux, flux_err = line.strip().split()
        if hit_type == 'exclusive':
            flux = float(flux)
            sex_intensities_given_read_name[read_name] = flux
    both_read_names = list(sex_intensities_given_read_name.keys())
    area_x = [area_given_read_name[read_name] for read_name in both_read_names]
    sextensity_y = [sex_intensities_given_read_name[read_name] for read_name in
                    both_read_names]
    r, pval = pearsonr(area_x, sextensity_y)
    plt.plot(area_x, sextensity_y, '.')
    plt.xlabel('Inferred Gaussian Area')
    plt.ylabel('Sextractor Flux')
    plt.savefig(fig_fpath)
    print 'Total time:', time.time() - start_time


def infer_all_nd2_point_areas(nd2_fpath, results_dir, fig_dir):
    """
    Infers all point areas in given nd2_fpath file.
    """
    def f(im_idx):
        infer_nd2_im_point_areas(nd2_fpath, im_idx, results_dir, fig_dir)
    from pathos.multiprocessing import ProcessingPool
    nd2 = nd2reader.Nd2(nd2_fpath)
    p = ProcessingPool(14)
    p.map(f, range(len(nd2)))
