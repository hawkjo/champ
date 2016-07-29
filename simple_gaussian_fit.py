import sys
import os
import argparse
import numpy as np
from scipy.optimize import minimize
import h5py
import hdf5_tools
import misc



def simple_gaussian(params, im_shape, im_center):
    """
    Calculate a simple gaussian

    `params` is a vector length 5 containing background intensity, amplitude, sigma, and
        mu_r, mu_c values for the gaussian 
    `im_shape` is the shape of the grid on which to calculate function values.
    """
    bg = params[0]
    A = params[1]
    inner_factor = -1.0 / (2 * params[2]**2)  # sigma = params[2] 
    mu_r = im_center[0] + params[3]
    mu_c = im_center[1] + params[4]

    r, c = np.arange(im_shape[0]), np.arange(im_shape[1])
    gr = np.exp(inner_factor * (r-mu_r)**2)
    gc = np.exp(inner_factor * (c-mu_c)**2)
    return bg + A * np.outer(gr, gc)


def fit_centered_gaussian(im, bg_bounds, A_bounds, max_del_mu, sigma_bounds):
    """
    Fit a simple gaussian approximately in the center of given image.
    """
    def gaussian_sq_err(params, im_shape, im_center, flattened_im):
        fit_g = simple_gaussian(params, im_shape, im_center)
        return sum((fit_val - im_val)**2 for fit_val, im_val in zip(fit_g.flatten(), flattened_im))

    flattened_im = im.flatten()
    im_shape = im.shape
    im_center = (int(im.shape[0]/2), int(im.shape[1]/2))

    bg0 = im.min()
    A0 = im.max() - bg0
    sigma0 = 1
    mu_r0, mu_c0 = 0, 0

    mu_r_bounds = (mu_r0 - max_del_mu, mu_r0 + max_del_mu)
    mu_c_bounds = (mu_c0 - max_del_mu, mu_c0 + max_del_mu)

    x0 = [bg0, A0, sigma0, mu_r0, mu_c0]
    bounds = [bg_bounds, A_bounds, sigma_bounds, mu_r_bounds, mu_c_bounds]
    args=(im_shape, im_center, flattened_im)

    res = minimize(gaussian_sq_err, x0, method='L-BFGS-B', args=args, bounds=bounds)
    if res.success:
        return res
    else:
        return False


def fit_gaussians_in_image(im, alignment_fpath, side_px, max_del_mu, sigma_bounds):
    """
    Fig centered gaussians on all aligned points in image
    """
    im_med = np.median(im.flatten())
    bg_bounds = (0.8 * im_med, 1.2 * im_med)
    A_bounds = (0.001 * im_med, None)
    rmax, cmax = im.shape[0] - side_px - 1, im.shape[1] - side_px - 1

    results = []
    for line in open(alignment_fpath):
        read_name, r, c = line.strip().split()
        r, c = map(misc.stoftoi, (r, c))
        if side_px <= r < rmax and side_px <= c < cmax:
            sub_im = im[r-side_px:r+side_px+1, c-side_px:c+side_px+1].astype(float)
            results.append((read_name, fit_centered_gaussian(sub_im, 
                                                             bg_bounds,
                                                             A_bounds,
                                                             max_del_mu, 
                                                             sigma_bounds)))
    return results


def fit_gaussians_in_h5_im(h5_fpath, results_dir, im_idx, side_px, max_del_mu, sigma_bounds):
    """
    Fit all gaussians in given h5 file image if aligned..
    """
    with h5py.File(h5_fpath) as f:
        channels = f.keys()
        for channel in channels:
            im_keys = f[channel].keys()
            im_keys.sort()
            im_key = im_keys[im_idx]
            bname = hdf5_tools.bname_given_channel_and_dset_name(channel, im_key)
            alignment_fname = '{}_all_read_rcs.txt'.format(bname)
            alignment_fpath = os.path.join(results_dir, alignment_fname)
            if not os.path.exists(alignment_fpath):
                sys.stdout.write('*')
                sys.stdout.flush()
                continue
            im = np.array(f[channel][im_key])

            results = fit_gaussians_in_image(im, 
                                             alignment_fpath, 
                                             side_px, 
                                             max_del_mu, 
                                             sigma_bounds)
            output_fpath = os.path.join(results_dir, '{}_gaussian_intensities.txt'.format(bname))
            with open(output_fpath, 'w') as out:
                for read_name, res in results:
                    if res:
                        out.write('\t'.join([read_name] + map(str, res.x)) + '\n')
                    else:
                        out.write('\t'.join([read_name] + ['-']*5) + '\n')
            sys.stdout.write('.')
            sys.stdout.flush()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('h5_fpath', help='Path to hdf5 file.')
    parser.add_argument('results_dir', help='Directory of alignment and output files.')
    parser.add_argument('im_idx', type=int)
    parser.add_argument('--side_px', action='store', dest='side_px', type=int, default=3,
                        help='Number of pixels to the side of center to include in fit.')
    parser.add_argument('--max_del_mu', action='store', dest='max_del_mu', type=float, default=1.0,
                        help='Maximum delta mu during fit.')
    parser.add_argument('--sigma_lb', action='store', dest='sigma_lb', type=float, default=1.0,
                        help='Lower bound of sigma')
    parser.add_argument('--sigma_ub', action='store', dest='sigma_ub', type=float, default=2.0,
                        help='Upper bound of sigma')

    args = parser.parse_args()
    fit_gaussians_in_h5_im(args.h5_fpath,
                           args.results_dir,
                           args.im_idx,
                           args.side_px,
                           args.max_del_mu,
                           (args.sigma_lb, args.sigma_ub),
                          )
