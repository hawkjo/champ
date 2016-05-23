from chimp.grid import GridImages
import functools
import logging
import multiprocessing
from multiprocessing import Pool
import os
import subprocess
import sys
import time
import h5py
import numpy as np


log = logging.getLogger(__name__)


class ImageFiles(object):
    def __init__(self, image_directory, filenames):
        self._image_directory = image_directory
        self._filenames = filenames

    def __len__(self):
        return len(self._filenames)

    # @property
    # def filenames(self):
    #     for f in self._filenames:
    #         yield os.path.join(self._image_directory, f)

    @property
    def directories(self):
        for f in self._filenames:
            yield os.path.join(self._image_directory, os.path.splitext(f)[0])


class XYZFile(object):
    def __init__(self, image):
        self._image = image

    def __str__(self):
        column_width = len(str(self._image.shape[0]))
        row_width = len(str(self._image.shape[1]))
        intensity_width = len(str(np.max(self._image)))
        line = "{column: >%s} {row: >%s} {intensity: >%s}" % (column_width, row_width, intensity_width)
        return "\n".join(line.format(row=row,
                                     column=column,
                                     intensity=intensity) for row, column, intensity in self._as_vector())

    def _as_vector(self):
        for (row, column), intensity in np.ndenumerate(self._image):
            yield row, column, intensity


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
        channel = str(channel).strip().replace(" ", "_")
        grid = GridImages(h5, channel)
        for n, image in enumerate(grid):
            xyz_file = XYZFile(image)
            xyz_path = "%s.xyz" % os.path.join(h5_base_name, image.index)
            with open(xyz_path, "w+") as f:
                f.write(str(xyz_file))
            fits_path = '%s.fits' % os.path.join(h5_base_name, image.index)
            subprocess.call(['fitsify', xyz_path, fits_path, '1', '2', '3'])
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
        # I think this is redundant since we created the HDF5 files from OME-TIFFs inside directories
        # that have this name already
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
