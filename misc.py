"""
A space for miscellaneous useful functions.
"""
import numpy as np
import re


def next_power_of_2(x):
    return 1 << (int(np.ceil(x)) - 1).bit_length()


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


def rcs_given_read_names(read_names):
    return np.array([map(int, name.split(':')[-2:]) for name in read_names])


def median_normalize(im):
    med = np.median(im)
    # Doing in place division by a float won't work because we have an int64 array
    # By casting to float with copy=False, we create a float view that allows
    # in place division without having to perform any copies
    im = im.astype('float', copy=False, casting='safe')
    im /= float(med)
    im -= 1.0
    return im


def miseq_version_given_read_names(read_names):
    tile_nums = set(int(name.split(':')[4]) for name in read_names)
    v2_tile_nums = set(range(1101, 1115) + range(2101, 2115))
    if tile_nums <= v2_tile_nums:
        return 'v2'
    else:
        return 'v3'


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


class AlignmentStats:
    def __init__(self, fpath):
        for line in open(fpath):
            name, value = line.strip().split(':')
            if '(' in name:
                name = name[:name.index('(')].strip()
            name = '_'.join(name.lower().split()).replace('-', '_')

            value = value.strip()
            if name == 'rc_offset':
                pattern = '\(([-.0-9]+),([-.0-9]+)\)'
                if value.count(',') > 1:
                    value = [np.array([float(m.group(1)), float(m.group(2))])
                             for m in re.finditer(pattern, value)]
                else:
                    m = re.match(pattern, value)
                    value = [np.array([float(m.group(1)), float(m.group(2))])]
            else:
                if ',' in value:
                    value = value.split(',')
                    if all([strisint(v) for v in value]):
                        value = map(int, value)
                    elif all([strisfloat(v) for v in value]):
                        value = map(float, value)
                else:
                    if strisint(value):
                        value = [int(value)]
                    elif strisfloat(value):
                        value = [float(value)]
                    else:
                        value = [value]
            setattr(self, name, value)

        self.numtiles = len(self.tile)


def pM_concentration_given_fpath(fpath):
    pattern = '[-_]([0-9_.]+)([pn]M)'
    m = re.search(pattern, fpath)
    assert m, fpath
    conc = float(m.group(1).replace('_', '.'))
    if m.group(2) == 'pM':
        return conc
    elif m.group(2) == 'nM':
        return conc * 1000
    else:
        raise ValueError('Can only handle pM and nM at the moment.')


def fold_radial_symmetry(x, with_max=False):
    """Takes square matrix and folds it by radial symmetry into summed entries, and optionally max."""
    # Matrix value types are summarized as follows, with 7x7 as example size:
    #
    #   d . . s . . d 
    #   . d . s . d . 
    #   . . d s d . . 
    #   s s s c s s s 
    #   . . d s d . . 
    #   . d . s . d . 
    #   d . . s . . d 
    #
    # Number of entries collapsed into one in this example via radial symmetry are as follows for
    # each value type:
    #
    #   c: 1
    #   s: 4
    #   d: 4
    #   .: 8
    #
    # This results in a reduction from 49 to 10 entries. With max gives 11 entries.

    slen = x.shape[0]
    assert x.shape == (slen, slen), x  # Check for 2d square matrix
    m = int(slen/2)
    folded = []
    
    # Center
    folded.append(x[m, m])
    
    # Sides and Diagonals
    for i in range(m):
        o = slen-i-1
        folded.append(x[i, m] + x[o, m] + x[m, i] + x[m, o])  # side 
        folded.append(x[i, i] + x[i, o] + x[o, i] + x[o, o])  # diagonal
        
    # Others
    for i in range(2, m+1):   # L_infty from c
        for j in range(1, i): # L_1 from s
            folded.append(sum([x[m-i, m-j],
                               x[m-i, m+j],
                               x[m+i, m-j],
                               x[m+i, m+j],
                               x[m-j, m-i],
                               x[m+j, m-i],
                               x[m-j, m+i],
                               x[m+j, m+i]]))

    # Max
    if with_max:
        folded.append(x.max())
    
    return folded


def unfold_radial_symmetry(folded, with_max=False):
    """Undo fold_radial_symmetry, where collapsed pixels are divided evenly to return them."""
    if isinstance(folded, np.ndarray):
        assert min(folded.shape) == 1, folded
        folded = list(folded.flatten())

    # The following determines the value of m as defined above
    m = (-1.5 + np.sqrt(1.5**2 - 4 * 0.5 * (1 - len(folded))))
    assert m == int(m), len(folded)
    m = int(m)
    slen = 2 * m + 1
    x = np.empty((slen, slen))

    # Center
    x[m, m] = folded.pop(0)

    # Sides and Diagonals
    for i in range(m):
        o = slen-i-1

        sval = folded.pop(0) / 4.0
        x[i, m] = sval
        x[o, m] = sval
        x[m, i] = sval
        x[m, o] = sval

        dval = folded.pop(0) / 4.0
        x[i, i] = dval
        x[i, o] = dval
        x[o, i] = dval
        x[o, o] = dval

    # Others
    for i in range(2, m+1):   # L_infty from c
        for j in range(1, i):  # L_1 from s
            val = folded.pop(0) / 8.0
            x[m-i, m-j] = val
            x[m-i, m+j] = val
            x[m+i, m-j] = val
            x[m+i, m+j] = val
            x[m-j, m-i] = val
            x[m+j, m-i] = val
            x[m-j, m+i] = val
            x[m+j, m+i] = val

    if with_max:
        assert len(folded) == 1, folded
        return x, folded.pop(0)
    else:
        assert len(folded) == 0, folded
        return x


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
