import logging
import multiprocessing
import os
import time
from champ.find_clusters import otsu, sextractor

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


def ensure_image_data_directory_exists(h5_filename):
    """
    Creates a directory based on the HDF5 filenames in order to store data derived from them.
    """
    new_directory = os.path.join(h5_filename)
    if not os.path.isdir(new_directory):
        os.mkdir(new_directory)


def run(image_directory, strategy_name):
    strategies = {'ostu': otsu,
                  'sextractor': sextractor}
    strategy = strategies[strategy_name]

    image_files = ImageFiles(image_directory,
                             [f for f in os.listdir(image_directory) if f.endswith('.h5')])
    for directory in image_files.directories:
        ensure_image_data_directory_exists(directory)
    # Try to use one core per file, but top out at the number of cores that the machine has.
    thread_count = min(len(image_files), multiprocessing.cpu_count() - 2)
    log.debug("Using %s threads for cluster finding." % thread_count)
    start = time.time()
    strategy.main(image_files, thread_count)
    log.info("Done with cluster location. Elapsed time: %s seconds" % round(time.time() - start, 0))
