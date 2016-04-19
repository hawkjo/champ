import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import nd2reader
import nd2tools
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


def sum_of_gaussians_given_separated_params(im_shape, As, sigs, mus, beta):
    """
    Calculate sum of gaussians.
    """
    params = [beta]
    for A, sig in zip(As, sigs):
        params.extend([sig, A])
    return sum_of_gaussians(mus, im_shape, params)
    

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


def taylor_fit_im(im,
                  stats_fpath,
                  rcs_fpath,
                  lambda_rho,
                  lambda_sigma,
                  lambda_mu,
                  lambda_beta,
                  results_dir,
                  fig_dir,
                  im_idx):
    start_time = time.time()

    # Check for good alignment
#    print 'Check alignment'
#    if (not os.path.isfile(stats_fpath)) or alignment_obviously_bad(stats_fpath, rcs_fpath):
#        return

    # Bin points
    print 'Binning points'
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
    amplitudes, sigmas, mus, areas = {}, {}, {}, {}
    stats = Counter()
    beta_0 = np.percentile(im, 5)
    for (rbin, cbin), bin_read_names in sorted(binned_read_names.items()):
        print rbin, cbin, len(bin_read_names),
        bin_points = binned_adjusted_points[(rbin, cbin)]
        bin_im = im[rbin[0]:rbin[1], cbin[0]:cbin[1]]
        sig_0s = np.ones((len(bin_read_names),))
        bin_As, bin_sigs, bin_mus, bin_beta = taylor_gaussian_fit(bin_im,
                                                                  bin_points,
                                                                  sig_0s,
                                                                  beta_0,
                                                                  lambda_rho,
                                                                  lambda_sigma,
                                                                  lambda_mu,
                                                                  lambda_beta)

        bin_areas = [2 * A * sig * np.pi for A, sig in zip(bin_As, bin_sigs)]

        assert len(bin_read_names) == len(bin_points) == len(bin_areas) > 0, (len(bin_read_names), len(bin_points), len(bin_areas))
        rmax = rbin[1] - rbin[0] - boundary
        cmax = cbin[1] - cbin[0] - boundary
        for read_name, point, A, sig, mu, area in zip(bin_read_names,
                                                      bin_points,
                                                      bin_As,
                                                      bin_sigs,
                                                      bin_mus,
                                                      bin_areas):
            if boundary <= point[0] < rmax and boundary <= point[1] < cmax:  # Discard boundary
                amplitudes[read_name] = A
                sigmas[read_name] = sig
                mus[read_name] = mu
                areas[read_name] = area

    for stat, count in stats.items():
        print 'Total %s: %d' % (stat, count)
    print 'Inferrence time:', time.time() - start_time

    # Write areas
    print 'Write areas'
    areas_fpath = os.path.join(results_dir, '%d_gaussian_areas.txt' % im_idx)
    with open(areas_fpath, 'w') as out:
        out.write('# read_name\tamplitude\tsigma\tmur\tmuc\tarea\n')
        out.write('\n'.join(['%s\t%g\t%g\t%g' % (rn, amplitudes[rn], sigmas[rn],
                                                 mus[rn][0], mus[rn][1], areas[rn])
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


def make_distance_matrices_function(im_shape):
    # Initialize complex matrix for distance calculations
    z = np.zeros(im_shape, dtype=complex)
    for i in range(z.shape[0]):
        z[i, :] += i
    for j in range(z.shape[1]):
        z[:, j] += complex(0, j)

    def distance_matrices_function(pt):
        rc_diffs = z - complex(pt[0], pt[1])
        l2_dists = abs(rc_diffs)
        r_diffs = rc_diffs.real
        c_diffs = rc_diffs.imag
        return l2_dists, r_diffs, c_diffs
    return distance_matrices_function


def taylor_gaussian_fit(im,
                        read_names,
                        points,
                        sig_0,
                        beta_0,
                        lambda_rho,
                        lambda_sigma,
                        lambda_mu,
                        lambda_beta):
    """
    Fits gaussians to image data via taylor approximation and dense linear least squares.

    Params
        im:             (2-D matrix) Image data. Must be reasonably small due to dense solver.
        points:         (k-by-2 matrix) Initial [row, column] guesses for gaussian centers.
        beta_0:         (float) Initial guess for background noise.
        lambda_rho:     (float) Regularization parameter for all rho = sqrt(A)
        lambda_sigma:   (float) Regularization parameter for all sigmas
        lambda_mu:      (float) Regularization parameter for all mus
        lambda_beta:    (float) Regularization parameter for beta
    """
    # Find initial guesses for rhos
    zero_read_names = set()
    A_0s, sig_0s, mu_0s = [], [], []
    for pt, read_name in zip(points, read_names):
        A = im[pt[0], pt[1]] - beta_0
        if A >= 0:
            A_0s.append(A)
            mu_0s.append(pt)
            sig_0s.append(sig_0)
        else:
            zero_read_names.add(read_name)
    mu_0s = np.array(mu_0s)

    rho_0s = [np.sqrt(A0) for A0 in A_0s]

    # im_minus_f0
    f_0 = sum_of_gaussians_given_separated_params(im.shape, A_0s, sig_0s, mu_0s, beta_0)

    im_minus_f0 = (im - f_0).flatten()
    
    # Derivatives matrix
    k = mu_0s.shape[0]  # Number of gaussians

    D = np.zeros((len(im_minus_f0), 4*k))  # Matrix of derivative coefficients for gaussians terms

    distance_matrices_function = make_distance_matrices_function(im.shape)

    for i, (rho, sig, pt) in enumerate(zip(rho_0s, sig_0s, mu_0s)):
        gauss = sum_of_gaussians_given_separated_params(im.shape, [1.0], [sig], mu_0s[i:i+1, :], 0)
        l2_dists, r_diffs, c_diffs = distance_matrices_function(pt)

        # df/d rho
        D[:, i] = gauss.flatten() * 2 * rho

        # df/d sigma
        D[:, i + k] = np.multiply(gauss, np.square(l2_dists)).flatten() * (rho**2 / sig**3)

        # df/d mur
        D[:, i + 2*k] = np.multiply(gauss, r_diffs).flatten() * (rho / sig)**2

        # df/d muc
        D[:, i + 3*k] = np.multiply(gauss, c_diffs).flatten() * (rho / sig)**2

    # Regularization
    n = D.shape[1] + 1
    assert n == 4 * k + 1, (n, k)

    Lbda_diag = np.ones((n,))
    Lbda_diag[:k] = lambda_rho
    Lbda_diag[k:2*k] = lambda_sigma
    Lbda_diag[2*k:4*k] = lambda_mu
    Lbda_diag[4*k] = lambda_beta

    # Build final matrix and vector. They have the form
    #
    #       [ D | 1 ]           [ im_minus_f0 ]
    #   A = [ ----- ]       b = [ ----------- ]
    #       [ Lbda  ],          [      0      ]
    #
    # Where Lbda is diagonal with regularization params

    A = np.r_[
            np.c_[D, np.ones((D.shape[0], 1))],
            np.diag(Lbda_diag)
    ]

    b = np.r_[
            im_minus_f0.reshape(D.shape[0], 1),
            np.zeros((n, 1))
    ]

    # Solve
    delta, residuals, rank, sing = np.linalg.lstsq(A, b)
    assert delta.shape == (4*k + 1, 1), (delta.shape, k)

    # Find final values
    rhos = np.array(rho_0s).reshape(len(rho_0s), 1) + delta[:k]
    As = np.square(rhos)
    sigs = sig_0s + delta[k:2*k]
    mus = mu_0s + np.c_[delta[2*k:3*k], delta[3*k:4*k]]
    beta = beta_0 + delta[4*k]

    # Add in zero-intensity points
    out_As, out_sigs, out_mus = [], [], []
    i = 0
    for read_name, pt in zip(read_names, points):
        if read_name in zero_read_names:
            out_As.append(0)
            out_sigs.append(sig_0)
            out_mus.append(pt)
        else:
            out_As.append(As[i, 0])
            out_sigs.append(sigs[i, 0])
            out_mus.append(mus[i])
            i += 1
    assert i == len(As), (i, len(As))
    out_mus = np.array(out_mus)

    return out_As, out_sigs, out_mus, beta, delta, residuals, rank, sing
