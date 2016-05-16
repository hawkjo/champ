import h5py
import re


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
