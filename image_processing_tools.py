"""
A space to create commonly used functions and classes for image processing.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def signal_hist_and_func(im, plot_title='', verbose=True, plot_curves=True):
    """Determines empirical background noise curve and returns (good signal)/(all signal)
    interpolation function."""

    hy, bin_edges = np.histogram(im, 1000, normed=True)
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
    if verbose:
        print 'Non-trivial cdf range: %f - %f' % (l_bnd, u_bnd)

    delta_x = np.mean(hx[1:]-hx[:-1])
    if verbose:
        print 'delta_x:', delta_x

    num_ext_points = 10
    extended_x = np.r_[[min(0, 1.1*hx[0])], hx, 1.1*hx[-1]]
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


