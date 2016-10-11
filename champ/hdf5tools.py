import re

import h5py


def parse_coordinates(position_key):
    m = re.search('\((\d+), (\d+)\)', position_key)
    column = int(m.group(1))
    row = int(m.group(2))
    return column, row


def get_image_key(row, column):
    return '(Major, minor) = ({}, {})'.format(column, row)


def load_channel_names(h5_path):
    with h5py.File(h5_path) as f:
        return f.keys()


def get_all_image_positions(h5_path):
    all_columns, all_rows = set(), set()
    with h5py.File(h5_path) as f:
        for channel in f.keys():
            for position_key in f[channel].keys():
                column, row = parse_coordinates(position_key)
                all_columns.add(column)
                all_rows.add(row)
    return all_columns, all_rows


def calculate_grid_dimensions(h5_path):
    num_columns, num_rows = map(len, get_all_image_positions(h5_path))
    return num_columns, num_rows
