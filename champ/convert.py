import tifffile
import os
import numpy as np
from champ.tiff import TifsPerConcentration, TifsPerFieldOfView
from champ import misc
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
            if not filename.endswith('.tif'):
                continue
            paths[directory].add(os.path.join(directory, filename))
    return paths


def main(paths, flipud, fliplr):
    image_adjustments = []
    if flipud:
        image_adjustments.append(lambda x: np.flipud(x))
    if fliplr:
        image_adjustments.append(lambda x: np.fliplr(x))
    with h5py.File('images.h5', 'a') as h5:
        for directory, tifs in paths.items():
            condition = misc.parse_concentration(directory)
            condition_group = h5.create_group(condition)
            tiff_stack = load_tiff_stack(list(tifs), image_adjustments)
            for t in tiff_stack:
                for channel, image in t:
                    if channel not in h5[condition]:
                        group = condition_group.create_group(channel)
                    else:
                        group = condition_group[channel]
                    dataset = group.create_dataset(t.dataset_name, image.shape, dtype=image.dtype)
                    dataset[...] = image
            log.debug("Done converting images for condition: %s" % condition)
