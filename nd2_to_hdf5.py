import tifffile
import h5py
import hdf5_tools
import nd2reader
import nd2tools
import sys
import os
import glob
import re
import numpy as np
from collections import defaultdict


def nd2_dir_to_hdf5(h5_fpath, nd2_fpath, channel_idxs):
    nd2 = nd2reader.Nd2(nd2_fpath)
    channel_names = [nd2.channels[idx] for idx in channel_idxs]
    assert len(channel_names) == len(set(channel_names)), 'Non-unique channel names in channels of interest.'

    coord_info, xs, ys, zs, pos_names, rows, cols = nd2tools.get_nd2_image_coord_info(nd2)
    rows.sort()
    if len(rows) > len(cols):
        Major_axis = 'rows'
    else:
        Major_axis = 'cols'
    
    print h5_fpath
    with h5py.File(h5_fpath, 'a') as f:
        for im_idx, image in enumerate(nd2):
            channel_idx = im_idx % len(nd2.channels)
            if channel_idx not in channel_idxs:
                continue
            sys.stdout.write('.')
            sys.stdout.flush()

            channel_name = nd2.channels[channel_idx]

            pos_name = nd2tools.convert_nd2_coordinates(nd2, outfmt='pos_name', im_idx=im_idx)
            row = rows.index(pos_name[0])
            col = int(pos_name[1:])
            if Major_axis == 'rows':
                Major_axis_pos, minor_axis_pos = row, col
            else:
                Major_axis_pos, minor_axis_pos = col, row

            out_im = image.data

            dset_name = hdf5_tools.dset_name_given_coords(Major_axis_pos, minor_axis_pos)

            if channel_name not in f:
                g = f.create_group(channel_name)
            else:
                g = f[channel_name]

            dset = g.create_dataset(dset_name, out_im.shape, dtype=out_im.dtype)
            dset[...] = out_im
    print
            

if __name__ == '__main__':
    usg_fmt = '{} <hdf5_fpath> <nd2_fpath> <channel_idx_1> [<channel_idx_2> ...]'.format(sys.argv[0])
    if len(sys.argv) < 4:
        sys.exit(usg_fmt)

    h5_fpath = sys.argv[1]
    nd2_fpath = sys.argv[2]

    channel_idxs = map(int, sys.argv[3:])

    nd2_dir_to_hdf5(h5_fpath, nd2_fpath, channel_idxs)
