import tifffile
import os
import numpy as np
from chimp.tiff import TifsPerConcentration, TifsPerFieldOfView
from collections import defaultdict
import h5py
import logging

log = logging.getLogger(__name__)


def load_tiff_stack(tifs, adjustments):
    # figure out if we have one tif per field of view and concentration,
    # or if each tif contains every image for every field of view in a single concentration
    # then put the files into the appropriate class
    tif = tifffile.TiffFile(tifs[0])
    if len(tif) > tif.micromanager_metadata['summary']['Positions']:
        # We have a single file that contains every image for an entire concentration
        return TifsPerConcentration(tifs, adjustments)
    # Each field of view is in its own tif
    return TifsPerFieldOfView(tifs, adjustments)


def get_all_tif_paths(root_directory):
    paths = defaultdict(set)
    for directory, subdirs, filenames in os.walk(root_directory):
        if not filenames:
            continue
        for filename in filenames:
            if not filename.endswith('.ome.tif'):
                continue
            paths[directory].add(os.path.join(directory, filename))
    return paths


def main(paths, flipud, fliplr):
    image_adjustments = []
    if flipud:
        image_adjustments.append(lambda x: np.flipud(x))
    if fliplr:
        image_adjustments.append(lambda x: np.fliplr(x))

    for directory, tifs in paths.items():
        hdf5_filename = directory + ".h5"
        with h5py.File(hdf5_filename, 'a') as h5:
            tiff_stack = load_tiff_stack(list(tifs), image_adjustments)
            for t in tiff_stack:
                for channel, image in t:
                    if channel not in h5:
                        group = h5.create_group(channel)
                    else:
                        group = h5[channel]
                    dataset = group.create_dataset(t.dataset_name, image.shape, dtype=image.dtype)
                    dataset[...] = image
        log.debug("Done with %s" % hdf5_filename)

# def tif_dir_to_hdf5(hdf5_file_path, tif_stack, flipud, fliplr):
#     """
#     Converts a set of OME-TIFF files to HDF5, the format that CHIMP uses.
#
#     Parameters
#     ----------
#     hdf5_file_path      Path of the file to create. Can be relative.
#     tif_file_paths      List of filenames. Can be relative.
#     flipud              Flip the raw image data across the horizontal axis.
#     fliplr              Flip the raw image data across the vertical axis.
#
#     """
#     # build up a function to perform the transformations needed to get the image in the orientation
#     # that CHIMP is expecting
#     image_adjustments = []
#     if flipud:
#         image_adjustments.append(lambda x: np.flipud(x))
#     if fliplr:
#         image_adjustments.append(lambda x: np.fliplr(x))
#     #
#     # with h5py.File(hdf5_file_path, 'a') as h5:
#     #     for file_path in tqdm.tqdm(tif_file_paths):
#     #         major_axis_position, minor_axis_position = tif_axes[file_path]
#     #
#     #         with tifffile.TiffFile(file_path) as tif:
#     #             summary = tif.micromanager_metadata['summary']
#     #
#     #             # Find channel names and assert unique
#     #             channel_names = [sanitize_name(name) for name in summary['ChNames']]
#     #             assert summary['Channels'] == len(channel_names) == len(set(channel_names)), channel_names
#     #
#     #             # channel_idxs map tif pages to channels
#     #             channel_indexes = tif.micromanager_metadata['index_map']['channel']
#     #
#     #             # Setup defaultdict
#     #             height, width = summary['Height'], summary['Width']
#     #             summed_images = defaultdict(lambda *x: np.zeros((height, width), dtype=np.int))
#     #
#     #             # Add images
#     #             for channel_index, page in zip(channel_indexes, tif.pages):
#     #                 image = page.asarray()
#     #                 for adjustment in image_adjustments:
#     #                     image = adjustment(image)
#     #                 summed_images[channel_index] += image
#
#     # Add images to hdf5
#     dataset_name = '(Major, minor) = ({}, {})'.format(major_axis_position, minor_axis_position)
#     for index, channel_name in enumerate(channel_names):
#         if channel_name not in h5:
#             group = h5.create_group(channel_name)
#         else:
#             group = h5[channel_name]
#
#         out_image = summed_images[index]
#         dataset = group.create_dataset(dataset_name, out_image.shape, dtype=out_image.dtype)
#         dataset[...] = out_image
#
#
