import re

import h5py


def coords_given_dset_name(dset_name):
    m = re.search('\((\d+), (\d+)\)', dset_name)
    Major_pos = int(m.group(1))
    minor_pos = int(m.group(2))
    return Major_pos, minor_pos


def Major_pos_given_dset_name(dset_name):
    return coords_given_dset_name(dset_name)[0]


def bname_given_channel_and_dset_name(channel, dset_name):
    Major_pos, minor_pos = coords_given_dset_name(dset_name)
    return 'Channel_{}_Pos_{}_{}'.format(channel, Major_pos, minor_pos)


def dset_name_given_coords(Major_pos, minor_pos):
    return '(Major, minor) = ({}, {})'.format(Major_pos, minor_pos)


def get_channel_names(h5_fpath):
    with h5py.File(h5_fpath) as f:
        return f.keys()


def get_all_Major_minor_pos(h5_fpath):
    all_Major_pos, all_minor_pos = set(), set()
    with h5py.File(h5_fpath) as f:
        for channel in f.keys():
            for dset_name in f[channel].keys():
                Major_pos, minor_pos = coords_given_dset_name(dset_name)
                all_Major_pos.add(Major_pos)
                all_minor_pos.add(minor_pos)
    return all_Major_pos, all_minor_pos


def get_nMajor_nminor_pos(h5_fpath):
    return map(len, get_all_Major_minor_pos(h5_fpath))
