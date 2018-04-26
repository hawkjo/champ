import tifffile
import os
import numpy as np
from champ.tiff import TifsPerConcentration, TifsPerFieldOfView, sanitize_name
from collections import defaultdict
import h5py
import logging
import time


log = logging.getLogger(__name__)


def load_channel_names(tifs):
    channels = set()
    for filename in tifs:
        tif = tifffile.TiffFile(filename)
        raw_channel_names = tif.micromanager_metadata['summary']['ChNames']
        if type(raw_channel_names) in (str, unicode):
            channel_names = [raw_channel_names]
        else:
            channel_names = raw_channel_names
        for channel in channel_names:
            channels.add(sanitize_name(channel))
    return tuple(channels)


def load_tiff_stack(tifs, adjustments, min_column, max_column):
    # figure out if we have one tif per field of view and concentration,
    # or if each tif contains every image for every field of view in a single concentration
    # then put the files into the appropriate class
    tif = tifffile.TiffFile(tifs[0])
    if len(tifs) == tif.micromanager_metadata['summary']['Positions']:
        # Each field of view is in its own tif. These may be tif stacks if multiple exposures were taken
        return TifsPerFieldOfView(tifs, adjustments, min_column, max_column)
    # We have a single file that contains every image for an entire concentration
    return TifsPerConcentration(tifs, adjustments, min_column, max_column)


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
                    if t.dataset_name not in group:
                        dataset = group.create_dataset(t.dataset_name, image.shape, dtype=image.dtype)
                    else:
                        dataset = group[t.dataset_name]
                    dataset[...] = image
        log.debug("Done with %s" % hdf5_filename)
