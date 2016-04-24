import tifffile
import h5py
import sys
import os
import glob
import re
import numpy as np


def tif_dir_to_hdf5(Major_axis_idx, channel, hdf5_fpath, tif_fpaths):
    Major_axis_group = Major_axis_idx + 1
    minor_axis_group = (1 - Major_axis_idx) + 1
    with h5py.File(hdf5_fpath, 'a') as f:
        g = f.create_group(channel)
        for fpath in tif_fpaths:
            sys.stdout.write('.')
            sys.stdout.flush()

            fname = os.path.split(fpath)[1]
            m = re.search('Pos_(\d+)_(\d+)\D', fname)
            Major_axis_pos = int(m.group(Major_axis_group))
            minor_axis_pos = int(m.group(minor_axis_group))

            im = tifffile.imread(fpath)
            if not (len(im.shape) == 3 and im.shape[1] == im.shape[2]):
                print '\nWarning: Skipping im_shape {} in {}'.format(im.shape, fpath)
                continue
            out_im = im.sum(axis=0)
            assert out_im.shape == im.shape[1:], out_im.shape

            dset_name = '(Major, minor) = ({}, {})'.format(Major_axis_pos, minor_axis_pos)
            dset = g.create_dataset(dset_name, out_im.shape, dtype=out_im.dtype)
            dset [...] = out_im
    print
            

if __name__ == '__main__':
    usg_fmt = '{} <Major_axis_idx(0,1)> <channel> <hdf5_fpath> <tif_fpaths>'.format(sys.argv[0])
    if len(sys.argv) < len(usg_fmt.split()):
        sys.exit(usg_fmt)

    Major_axis_idx = int(sys.argv[1])
    assert Major_axis_idx in [0, 1], 'Major axis must be in (0, 1)'

    channel = sys.argv[2]

    hdf5_fpath = sys.argv[3]

    tif_fpaths = sys.argv[4:]

    tif_dir_to_hdf5(Major_axis_idx, channel, hdf5_fpath, tif_fpaths)
