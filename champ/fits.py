from champ.grid import GridImages
import functools
import logging
import multiprocessing
from multiprocessing import Pool
import os
import subprocess
import sys
import time
import h5py
from astropy.io import fits
import numpy as np


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


class SEConfig(object):
    def __enter__(self):
        self._create_config_files()

    def __exit__(self, *args):
        self._delete_config_files()

    def _delete_config_files(self):
        for filename in ('default.sex', 'spot.param', 'default.conv'):
            try:
                os.unlink(filename)
            except OSError:
                pass

    def _create_config_files(self):
        default_text = """DETECT_THRESH 2
DEBLEND_NTHRESH 64
DEBLEND_MINCONT 0.00005
"""
        with open('default.sex', 'w+') as f:
            f.write(default_text)

        spot_text = """X_IMAGE
Y_IMAGE
FLUX_AUTO
FLUXERR_AUTO
FLAGS
A_IMAGE
B_IMAGE
THETA_IMAGE
"""
        with open('spot.param', 'w+') as f:
            f.write(spot_text)

        convolution_text = """CONV NORM
1 2 1
2 4 2
1 2 1
"""
        with open('default.conv', 'w+') as f:
            f.write(convolution_text)


def get_base_file_names(h5_filename):
    return ["%s" % os.path.join(h5_filename, os.path.splitext(filename)[0])
            for filename in os.listdir(h5_filename) if filename.endswith(".xyz")]


def source_extract(base_file):
    command = '/usr/bin/sextractor {base_file}.fits -PARAMETERS_NAME spot.param -CATALOG_NAME {base_file}.cat -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME {base_file}.model'
    # Don't print any output
    with open('/dev/null', 'w') as devnull:
        command = command.format(base_file=base_file).split(' ')
        subprocess.call(command, stdout=devnull, stderr=devnull)


def create_fits_files(h5_base_name):
    h5_filename = h5_base_name + ".h5"
    log.info("Creating fits files for %s" % h5_filename)
    h5 = h5py.File(h5_filename)
    for channel in h5.keys():
        grid = GridImages(h5, channel)
        for image in grid:
            fits_path = '%s.fits' % os.path.join(h5_base_name, image.index)
            # Source Extractor can handle at most 32-bit values, so we have to cast down from our 64-bit images or
            # else it will throw an error. We clip to ensure there's no overflow, although this seems improbable
            # given that most cameras are 16 bit
            clipped_image = np.clip(image, 0, 2**32-1).astype(np.uint32)
            hdu = fits.PrimaryHDU(clipped_image)
            hdu.writeto(fits_path, clobber=True)
    log.info("Done creating fits files for %s" % h5_base_name)


def ensure_image_data_directory_exists(h5_filename):
    """
    Creates a directory based on the HDF5 filenames in order to store data derived from them.

    """
    new_directory = os.path.join(h5_filename)
    if not os.path.isdir(new_directory):
        os.mkdir(new_directory)


def main(image_directory):
    image_files = ImageFiles(image_directory,
                             [f for f in os.listdir(image_directory) if f.endswith('.h5')])
    for directory in image_files.directories:
        ensure_image_data_directory_exists(directory)
    # Try to use one core per file, but top out at the number of cores that the machine has.
    thread_count = min(len(image_files), multiprocessing.cpu_count())
    log.debug("Using %s threads for source extraction" % thread_count)
    # Assign each HDF5 file to a thread, which converts it to a "fits" file
    worker_pool = Pool(processes=thread_count)

    log.info("Starting fits file conversions.")
    # The multiprocessing thing only takes an iterable with no arguments, so we use a partial function to pass
    # the directory where the files should be written
    fits_func = functools.partial(create_fits_files)
    # KeyboardInterrupt won't behave as expected while multiprocessing unless you specify a timeout.
    # We don't want one really, so we just use the largest possible integer instead
    start = time.time()
    worker_pool.map_async(fits_func, image_files.directories).get(timeout=sys.maxint)
    log.info("Done with fits file conversions. Elapsed time: %s seconds" % round(time.time() - start, 0))

    # Now run source extractor to find the coordinates of points
    with SEConfig():
        log.info("Starting Source Extractor...")
        start = time.time()
        # Set up a worker for each HDF5 file like before
        worker_pool = Pool(thread_count)
        base_files = [base_file for h5_filename in image_files.directories
                      for base_file in get_base_file_names(h5_filename)]
        worker_pool.map_async(source_extract, base_files).get(timeout=sys.maxint)
        log.info("Done with Source Extractor! Took %s seconds" % round(time.time() - start, 0))
