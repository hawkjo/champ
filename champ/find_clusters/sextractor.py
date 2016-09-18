from astropy.io import fits
from champ.grid import GridImages
import logging
import numpy as np
import os
from multiprocessing import Pool
import subprocess
import sys
import time
import h5py
import functools

log = logging.getLogger(__name__)


def main(image_files, thread_count):
    fits_func = functools.partial(create_fits_files)
    # KeyboardInterrupt won't behave as expected while multiprocessing unless you specify a timeout.
    # We don't want one really, so we just use the largest possible integer instead
    start = time.time()
    worker_pool = Pool(thread_count)
    worker_pool.map_async(fits_func, image_files.directories).get(timeout=sys.maxint)
    log.info("Done with FITS files. Time: %s" % (time.time() - start))
    with SEConfig():
        log.info("Starting Source Extractor...")
        start = time.time()
        # Set up a worker for each HDF5 file like before
        worker_pool = Pool(thread_count)
        base_files = [base_file for h5_filename in image_files.directories
                      for base_file in get_base_file_names(h5_filename)]
        worker_pool.map_async(source_extract, base_files).get(timeout=sys.maxint)
        log.info("Done with Source Extractor! Took %s seconds" % round(time.time() - start, 0))


def get_base_file_names(h5_filename):
    return ["%s" % os.path.join(h5_filename, os.path.splitext(filename)[0])
            for filename in os.listdir(h5_filename) if filename.endswith(".fits")]


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


def source_extract(base_file):
    command = '/usr/bin/sextractor {base_file}.fits -PARAMETERS_NAME spot.param -CATALOG_NAME {base_file}.cat -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME {base_file}.model'
    # Don't print any output
    with open('/dev/null', 'w') as devnull:
        command = command.format(base_file=base_file).split(' ')
        subprocess.call(command, stdout=devnull, stderr=devnull)


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
