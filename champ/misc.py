"""
A space for miscellaneous useful functions.
"""
import re
import numpy as np
from sklearn.neighbors import KernelDensity
from scipy.optimize import minimize


def next_power_of_2(x):
    return 1 << (int(np.ceil(x))-1).bit_length()


def max_2d_idx(a):
    return np.unravel_index(a.argmax(), a.shape)


def pad_to_size(M, size):
    assert len(size) == 2, 'Row and column sizes needed.'
    left_to_pad = size - np.array(M.shape)
    return np.pad(M, ((0, left_to_pad[0]), (0, left_to_pad[1])), mode='constant')


def right_rotation_matrix(angle, degrees=True):
    if degrees:
        angle *= np.pi / 180.0
    sina = np.sin(angle)
    cosa = np.cos(angle)
    return np.array([[cosa, sina],
                     [-sina, cosa]])


def strisfloat(x):
    try:
        a = float(x)
    except ValueError:
        return False
    else:
        return True


def strisint(x):
    try:
        a = float(x)
        b = int(x)
    except ValueError:
        return False
    else:
        return a == b


def stoftoi(s):
    return int(round(float(s)))


def parse_concentration(filename):
    pattern = '[-_]([0-9_.]+)([pn]M)'
    m = re.search(pattern, filename)
    if m is None:
        raise ValueError("The concentration cannot be parsed from the filename: %s" % filename)
    conc = float(m.group(1).replace('_', '.'))
    if m.group(2) == 'pM':
        return conc
    elif m.group(2) == 'nM':
        return conc * 1000
    else:
        raise ValueError('Can only handle pM and nM at the moment.')


def read_names_and_points_given_rcs_fpath(rcs_fpath):
    """
    Return the read names and (r, c) point locations of points in implied image.
    """
    read_names, points = [], []
    for line in open(rcs_fpath):
        var = line.strip().split()
        read_names.append(var[0])
        points.append(map(float, var[1:]))
    return read_names, np.array(points)


def list_if_scalar(x, list_len):
    try:
        float(x)
        return [x]*list_len
    except:
        return


def get_mode(vals):
    h = 1.06 * np.std(vals) * len(vals)**(-1.0/5.0)
    kdf = KernelDensity(bandwidth=h)
    kdf.fit(np.array(vals).reshape(len(vals), 1))

    def neg_kdf(x):
        return -kdf.score(x)
    res = minimize(neg_kdf, x0=np.median(vals), method='Nelder-Mead')
    assert res.success, res
    return float(res.x)
