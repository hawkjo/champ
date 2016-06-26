import tifffile
import h5py
import os
import re
import numpy as np
from collections import defaultdict
import tqdm

whitespace_regex = re.compile('[\s]+')
special_chars_regex = re.compile('[\W]+')
name_regex = re.compile('Pos_(\d+)_(\d+)\D')


def tif_dir_to_hdf5(hdf5_file_path, tif_file_paths, flipud, fliplr):
    """
    Converts a set of OME-TIFF files to HDF5, the format that CHIMP uses.

    Parameters
    ----------
    hdf5_file_path      Path of the file to create. Can be relative.
    tif_file_paths      List of filenames. Can be relative.
    flipud              Flip the raw image data across the horizontal axis.
    fliplr              Flip the raw image data across the vertical axis.

    """
    # build up a function to perform the transformations needed to get the image in the orientation
    # that CHIMP is expecting
    image_adjustments = []
    if flipud:
        image_adjustments.append(lambda x: np.flipud(x))
    if fliplr:
        image_adjustments.append(lambda x: np.fliplr(x))

    tif_axes = build_tif_axes(tif_file_paths)
    with h5py.File(hdf5_file_path, 'a') as h5:
        for file_path in tqdm.tqdm(tif_file_paths):
            major_axis_position, minor_axis_position = tif_axes[file_path]

            with tifffile.TiffFile(file_path) as tif:
                summary = tif.micromanager_metadata['summary']

                # Find channel names and assert unique
                channel_names = [sanitize_name(name) for name in summary['ChNames']]
                assert summary['Channels'] == len(channel_names) == len(set(channel_names)), channel_names

                # channel_idxs map tif pages to channels
                channel_indexes = tif.micromanager_metadata['index_map']['channel']

                # Setup defaultdict
                height, width = summary['Height'], summary['Width']
                summed_images = defaultdict(lambda *x: np.zeros((height, width), dtype=np.int))

                # Add images
                for channel_index, page in zip(channel_indexes, tif.pages):
                    image = page.asarray()
                    for adjustment in image_adjustments:
                        image = adjustment(image)
                    summed_images[channel_index] += image

                # Add images to hdf5
                dataset_name = '(Major, minor) = ({}, {})'.format(major_axis_position, minor_axis_position)
                for index, channel_name in enumerate(channel_names):
                    if channel_name not in h5:
                        group = h5.create_group(channel_name)
                    else:
                        group = h5[channel_name]

                    out_image = summed_images[index]
                    dataset = group.create_dataset(dataset_name, out_image.shape, dtype=out_image.dtype)
                    dataset[...] = out_image


def sanitize_name(name):
    return special_chars_regex.sub('', whitespace_regex.sub('_', name))


def build_tif_axes(tif_file_paths):
    # we need to know the major and minor axis of each tif file, since we assume everywhere else in the codebase that
    # the Illumina tiles are arranged from left-to-right. However, depending on how you acquired your images, it could
    # be long up and down instead. We build up a dictionary of each tif file and its axes, guaranteeing
    # that the major axis comes first.
    tif_axes = {}
    best_first = 0
    best_second = 0
    for file_path in tif_file_paths:
        filename = os.path.split(file_path)[1]
        axis_positions = name_regex.search(filename)
        first = int(axis_positions.group(1))
        second = int(axis_positions.group(2))
        best_first = max(first, best_first)
        best_second = max(second, best_second)
        tif_axes[file_path] = (first, second)
    if best_second > best_first:
        # the second thing is the major axis, so we need to invert them
        return {file_path: (second, first) for file_path, (first, second) in tif_axes.items()}
    # no need to invert, just return the values we already have
    return tif_axes
