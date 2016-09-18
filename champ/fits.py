from champ.grid import GridImages
import logging
import multiprocessing
from multiprocessing import Pool
import os
import sys
import time
import h5py
from skimage.filters import threshold_otsu
from scipy import ndimage

log = logging.getLogger(__name__)


class ImageFiles(object):
    def __init__(self, image_directory, filenames):
        self._image_directory = image_directory
        self._filenames = filenames

    def __len__(self):
        return len(self._filenames)

    @property
    def directories(self):
        for f in self._filenames:
            yield os.path.join(self._image_directory, os.path.splitext(f)[0])


def get_base_file_names(h5_filename):
    return ["%s" % os.path.join(h5_filename, os.path.splitext(filename)[0])
            for filename in os.listdir(h5_filename) if filename.endswith(".xyz")]


def ensure_image_data_directory_exists(h5_filename):
    """
    Creates a directory based on the HDF5 filenames in order to store data derived from them.
    """
    new_directory = os.path.join(h5_filename)
    if not os.path.isdir(new_directory):
        os.mkdir(new_directory)


def find_clusters(h5_base_name):
    h5_filename = h5_base_name + ".h5"
    log.info("Finding clusters for %s" % h5_filename)
    h5 = h5py.File(h5_filename)
    for channel in h5.keys():
        grid = GridImages(h5, channel)
        for image in grid:
            out_filepath = os.path.join(h5_base_name, image.index + '.cat')
            threshold = threshold_otsu(image)
            mask_pixels = (image > threshold)
            mask = ndimage.binary_closing(ndimage.binary_opening(mask_pixels))
            label_image, num_labels = ndimage.label(mask)
            log.debug("Found %d clusters in %s/%s" % (num_labels, h5_base_name, image.index))
            center_of_masses = ndimage.center_of_mass(image, label_image, range(num_labels + 1))
            write_cluster_locations(center_of_masses, out_filepath)


def write_cluster_locations(locations, out_filepath):
    with open(out_filepath, 'w') as out:
        out.write('\n'.join("%s\t%s" % (r, c) for c, r in locations))


def main(image_directory):
    image_files = ImageFiles(image_directory,
                             [f for f in os.listdir(image_directory) if f.endswith('.h5')])
    for directory in image_files.directories:
        ensure_image_data_directory_exists(directory)
    # Try to use one core per file, but top out at the number of cores that the machine has.
    thread_count = min(len(image_files), multiprocessing.cpu_count() - 2)
    log.debug("Using %s threads for source extraction" % thread_count)
    # Assign each HDF5 file to a thread, which converts it to a "fits" file
    worker_pool = Pool(processes=thread_count)
    start = time.time()
    worker_pool.map_async(find_clusters, image_files.directories).get(timeout=sys.maxint)
    log.info("Done with cluster location. Elapsed time: %s seconds" % round(time.time() - start, 0))