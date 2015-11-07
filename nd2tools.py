import sys
import os
import copy
import nd2reader
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from misc import median_normalize


def get_nd2_image_coord_info(nd2):
    try:
        coord_info = nd2.metadata['ImageMetadata']['SLxExperiment']['uLoopPars']['Points']['']
    except:
        try:
            coord_info = nd2.metadata['ImageMetadata']['SLxExperiment']['ppNextLevelEx']['']['uLoopPars']['Points']['']
        except:
            raise
    xs = [pt['dPosX'] for pt in coord_info]
    ys = [pt['dPosY'] for pt in coord_info]
    zs = [pt['dPosZ'] for pt in coord_info]
    pos_names = [pt['dPosName'] for pt in coord_info]
    rows = list(sorted(set(name[0] for name in pos_names)))
    cols = list(sorted(set(name[1:] for name in pos_names)))
    return coord_info, xs, ys, zs, pos_names, rows, cols


def nrows_and_ncols(nd2):
    coord_info, xs, ys, zs, pos_names, rows, cols = get_nd2_image_coord_info(nd2)
    return len(rows), len(cols)


def get_ims_in_run(nd2):
    coord_info, xs, ys, zs, pos_names, rows, cols = get_nd2_image_coord_info(nd2)
    return len(pos_names) * len(nd2.channels)


def plot_nd2_grid(nd2, edge, channel, idx_start=None, idx_end=None, suptitle=''):
    assert channel < len(nd2.channels), (channel, len(nd2.channels))
    coord_info, xs, ys, zs, pos_names, rows, cols = get_nd2_image_coord_info(nd2)
    nrows = len(rows)
    ncols = len(cols)
    
    gs_kwargs = dict(wspace=0.01, hspace=0.01)

    fig, axs = plt.subplots(nrows, ncols, 
                            figsize=(ncols * edge, nrows * edge), 
                            gridspec_kw=gs_kwargs)
    mn = float('inf')
    mx = float('-inf')
    
    if idx_start is None:
        idx_start = 0
    if idx_end is None:
        idx_end = len(nd2)
    
    # Force right channel
    if (idx_start - channel) % len(nd2.channels) != 0:
        idx_start = int(idx_start / len(nd2.channels)) * len(nd2.channels) + channel
        
    for i in range(idx_start, idx_end, len(nd2.channels)):
        im = median_normalize(nd2[i].data)
        im_pos_idx = int(i / len(nd2.channels)) % (nrows * ncols)
        pos_name = pos_names[im_pos_idx]
        r = nrows - 1 - rows.index(pos_name[0])
        c = ncols - 1 - cols.index(pos_name[1:])
        mn = min(mn, im.min())
        mx = max(mx, im.max())
        axs[r, c].matshow(im, vmin=0, vmax=2)#, cmap=plt.get_cmap('Blues'))
        axs[r, c].set_xticks([])
        axs[r, c].set_yticks([])
        axs[r, c].text(im.shape[1]/2,
                       im.shape[0]/2,
                       pos_name,
                       color='white',
                       fontsize=24,
                       ha='center',
                       va='center')

    fig.suptitle(suptitle, fontsize=30)
    
    return mn, mx, fig


def stitch_nd2(nd2, channel, row_start=None, row_end=None, col_start=None, col_end=None, fast=True):
    coord_info, xs, ys, zs, pos_names, rows, cols = get_nd2_image_coord_info(nd2)
    used_rows = list(sorted(set(ys)))[row_start:row_end]
    used_cols = list(sorted(set(xs)))[col_start:col_end]

    im_shape = nd2[0].data.shape

    def im_idx_given_row_col(row, col):
        return zip(ys, xs).index((row, col)) * len(nd2.channels) + channel
        #return pos_names.index(row + col) * len(nd2.channels) + channel

    def best_alignment_offset(im1_idx, im2_idx, max_overlap=100, min_overlap=3, max_lateral=60, direction='lr'):
        assert direction in ['lr', 'ud'], direction
        im1 = median_normalize(nd2[im1_idx].data)
        im2 = median_normalize(nd2[im2_idx].data)
        assert im1.shape == im2.shape == im_shape, (im1.shape, im2.shape)
        max_score, best_overlap, best_lateral = float('-inf'), None, None
        for ovr in range(min_overlap, max_overlap+1):
            for lat in range(-max_lateral, max_lateral):
                if direction == 'lr':
                    if lat < 0:
                        score_mat = np.multiply(im1[:im_shape[0]+lat, im_shape[1]-ovr:],
                                                im2[-lat:, :ovr])
                    else:
                        score_mat = np.multiply(im1[lat:, im_shape[1]-ovr:],
                                                im2[:im_shape[0]-lat, :ovr])
                else:
                    if lat < 0:
                        score_mat = np.multiply(im1[im_shape[0]-ovr:, :im_shape[1]+lat],
                                                im2[:ovr, -lat:])
                    else:
                        score_mat = np.multiply(im1[im_shape[0]-ovr:, lat:],
                                                im2[:ovr, :im_shape[1]-lat])
                score = float(score_mat.sum()) / (score_mat.shape[0] * score_mat.shape[1])
                if score > max_score:
                    max_score, best_overlap, best_lateral = score, ovr, lat
        if direction == 'lr':
            offset = np.array([best_lateral, im_shape[1]-best_overlap])
        else:
            offset = np.array([im_shape[0]-best_overlap, best_lateral])
        return offset

    def best_alignment_pos(nascent_im,
                           im_idx,
                           curr_pos,
                           direction,
                           max_overlap=100,
                           min_overlap=5):
        im = median_normalize(nd2[im_idx].data)
        assert im.shape == im_shape, im.shape
        assert (len(direction) == 2
                and len(set(direction) & set('lr')) == 1
                and len(set(direction) & set('ud')) == 1), direction
        d1, d2 = direction
        half_max = int(max_overlap / 2)
        offset = np.array([0, 0])
        if d1 == 'l':
            offset[1] = -im_shape[1] + min_overlap
        elif d1 == 'r':
            offset[1] = im_shape[1] - max_overlap
        elif d1 == 'u':
            offset[0] = -im_shape[0] + min_overlap
        else:
            offset[0] = im_shape[0] - max_overlap

        if d2 in 'lr':
            offset[1] = -half_max
        else:
            offset[0] = -half_max
        assert 0 not in offset, (direction, offset)

        start_pos = curr_pos + offset
        max_score, best_pos = float('-inf'), None
        for r_off in range(min_overlap, max_overlap+1):
            for c_off in range(min_overlap, max_overlap+1):
                pos = start_pos + np.array([r_off, c_off])
                r, c = pos
                score_mat = np.multiply(im, nascent_im[r:r+im_shape[0], c:c+im_shape[1]])
                score = float(score_mat.sum()) / (np.count_nonzero(score_mat) + 1)
                if score > max_score:
                    max_score, best_pos = score, copy.deepcopy(pos)
        r, c = best_pos
        nascent_im[r:r+im_shape[0], c:c+im_shape[1]] = im
        return best_pos

    if fast:
        positions = {}
        for i, row in enumerate(used_rows):
            for j, col in enumerate(used_cols):
                sys.stdout.write('.')
                sys.stdout.flush()
                im2_idx = im_idx_given_row_col(row, col)
                if i == j == 0:
                    curr_pos = np.array([0, 0])
                elif j == 0:
                    im1_idx = im_idx_given_row_col(used_rows[i-1], col)
                    curr_pos = copy.deepcopy(positions[(used_rows[i-1], col)])  # Reset to beginning of previous row
                    curr_pos += best_alignment_offset(im1_idx, im2_idx, direction='ud')
                else:
                    im1_idx = im_idx_given_row_col(row, used_cols[j-1])
                    curr_pos += best_alignment_offset(im1_idx, im2_idx, direction='lr')
                positions[(row, col)] = copy.deepcopy(curr_pos)
    else:
        positions = {}
        nascent_im = np.zeros((im_shape[0] * (2 + len(used_rows)), im_shape[1] * (2 + len(used_cols))))
        for i, row in enumerate(used_rows):
            for j, col in enumerate(used_cols):
                sys.stdout.write('.')
                sys.stdout.flush()
                im_idx = im_idx_given_row_col(row, col)
                if i == j == 0:
                    curr_pos = np.array(im_shape)
                elif j == 0:
                    curr_pos = copy.deepcopy(positions[(used_rows[i-1], col)])  # Reset to beginning of previous row
                    curr_pos = best_alignment_pos(nascent_im, im_idx, curr_pos, direction='dr')
                else:
                    curr_pos = best_alignment_pos(nascent_im, im_idx, curr_pos, direction='rd')
                r, c = curr_pos
                positions[(row, col)] = copy.deepcopy(curr_pos)
    

    minr = min(pos[0] for pos in positions.values())
    minc = min(pos[1] for pos in positions.values())
    offset = np.array([minr, minc])
    for k in positions.keys():
        positions[k] = positions[k] - offset
    maxr = max(pos[0] for pos in positions.values()) + im_shape[0]
    maxc = max(pos[1] for pos in positions.values()) + im_shape[1]

    stch_im = np.empty((maxr, maxc))
    diff_im = np.zeros((maxr, maxc))
    for i, row in enumerate(used_rows):
        for j, col in enumerate(used_cols):
            im_idx = im_idx_given_row_col(row, col)
            im = median_normalize(nd2[im_idx].data)
            r, c = positions[(row, col)]
            stch_im[r:r+im.shape[0], c:c+im.shape[1]] = im
            sign = (-1)**(i+j)
            diff_im[r:r+im.shape[0], c:c+im.shape[1]] += sign*im

    return stch_im, diff_im, positions


def convert_nd2_coordinates(nd2, outfmt, **kwargs):
    systems = ['pos_name', 'im_idx', 'pos_idx', 'pos_coords']
    assert len(set(systems) & set(kwargs.keys())) == 1, 'Must give exactly one input in %s' % systems
    assert outfmt in systems, outfmt
    if 'run' not in kwargs:
        run = 0
    ims_in_run = get_ims_in_run(nd2)

    coord_info, xs, ys, zs, pos_names, rows, cols = get_nd2_image_coord_info(nd2)
    if 'im_idx' in kwargs:
        im_idx = kwargs['im_idx']
        pos_idx = int((im_idx % ims_in_run) / len(nd2.channels))
    elif 'pos_name' in kwargs:
        pos_name = kwargs['pos_name']
        pos_idx = pos_names.index(pos_name)
    elif 'pos_coords' in kwargs:
        pos_idx = zip(ys, xs).index(kwargs['pos_coords'])
    else:
        pos_idx = kwargs['pos_idx']
        
    if outfmt == 'pos_idx':
        return pos_idx
    elif outfmt == 'pos_name':
        return pos_names[pos_idx]
    elif outfmt == 'pos_coords':
        return zip(ys, xs)[pos_idx]
    elif outfmt == 'im_idx':
        assert kwargs['channel'] < len(nd2.channels), (kwargs['channel'], nd2.channels)
        return run * ims_in_run + pos_idx * len(nd2.channels) + kwargs['channel']
    else:
        raise ValueError('"%s" is not a recognized outfmt.' % outfmt)


def bname_given_nd2(nd2):
    return bname_given_fpath(nd2._filename)


def bname_given_fpath(fpath):
    return os.path.splitext(os.path.basename(fpath))[0]
