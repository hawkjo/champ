import tifffile
import os
import numpy as np
from champ.tiff import TifsPerConcentration, TifsPerFieldOfView, sanitize_name
from collections import defaultdict
import h5py
import logging


log = logging.getLogger(__name__)


def load_channel_names(tifs):
    channels = set()
    for filename in tifs:
        tif = tifffile.TiffFile(filename)
        for channel in tif.micromanager_metadata['summary']['ChNames']:
            channels.add(sanitize_name(channel))
    return tuple(channels)


def load_tiff_stack(tifs, adjustments, min_column, max_column):
    # figure out if we have one tif per field of view and concentration,
    # or if each tif contains every image for every field of view in a single concentration
    # then put the files into the appropriate class
    tif = tifffile.TiffFile(tifs[0])
    if len(tif) > tif.micromanager_metadata['summary']['Positions']:
        # We have a single file that contains every image for an entire concentration
        print("** tifs:", tifs)
        return TifsPerConcentration(tifs, adjustments, min_column, max_column)
    # Each field of view is in its own tif
    return TifsPerFieldOfView(tifs, adjustments, min_column, max_column)


def get_all_tif_paths(root_directory):
    paths = defaultdict(set)
    for directory, subdirs, filenames in os.walk(root_directory):
        if not filenames:
            continue
        for filename in filenames:
            if not filename.endswith('.tif'):
                continue
            paths[directory].add(os.path.join(directory, filename))
    return paths


def main(paths, flipud, fliplr, min_column, max_column):
    image_adjustments = []
    if flipud:
        image_adjustments.append(lambda x: np.flipud(x))
    if fliplr:
        image_adjustments.append(lambda x: np.fliplr(x))

    for directory, tifs in paths.items():
        hdf5_filename = directory + ".h5"
        if os.path.exists(hdf5_filename):
            log.warn("HDF5 file already exists, skipping creation: %s" % hdf5_filename)
            continue
        with h5py.File(hdf5_filename, 'a') as h5:
            tiff_stack = load_tiff_stack(list(tifs), image_adjustments, min_column, max_column)
            for t in tiff_stack:
                for channel, image in t:
                    if channel not in h5:
                        group = h5.create_group(channel)
                    else:
                        group = h5[channel]
                    dataset = group.create_dataset(t.dataset_name, image.shape, dtype=image.dtype)
                    dataset[...] = image
        log.debug("Done with %s" % hdf5_filename)
