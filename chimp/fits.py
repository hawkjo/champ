from chimp import files
from chimp.model.xyz import XYZFile
import functools
import logging
import multiprocessing
from multiprocessing import Pool
from nd2reader import Nd2
import os
import subprocess
import sys
import time


log = logging.getLogger(__name__)


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


def base_files(nd2_filename):
    return ["%s" % os.path.join(nd2_filename, os.path.splitext(filename)[0])
            for filename in os.listdir(nd2_filename) if filename.endswith(".xyz")]


def source_extract(base_file):
    command = '/usr/bin/sextractor {base_file}.fits -PARAMETERS_NAME spot.param -CATALOG_NAME {base_file}.cat -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME {base_file}.model'
    # Don't print any output
    with open('/dev/null', 'w') as devnull:
        command = command.format(base_file=base_file).split(' ')
        subprocess.call(command, stdout=devnull, stderr=devnull)


def create_fits_files(nd2_filename):
    log.info("Creating fits files for %s..." % nd2_filename)
    nd2 = Nd2(nd2_filename + ".nd2")
    for n, image in enumerate(nd2):
        xyz_file = XYZFile(image)
        xyz_path = "%s.xyz" % os.path.join(nd2_filename, str(n))
        with open(xyz_path, "w+") as f:
            f.write(str(xyz_file))
        fits_path = '%s.fits' % os.path.join(nd2_filename, str(n))
        subprocess.call(['fitsify', xyz_path, fits_path, '1', '2', '3'])
    log.info("Done creating fits files for %s" % nd2_filename)


def main():
    image_files = files.load_image_files()
    for directory in image_files.directories:
        files.ensure_image_data_directory_exists(directory)
    # Try to use one core per file, but top out at the number of cores that the machine has.
    thread_count = min(len(image_files), multiprocessing.cpu_count())
    log.debug("Using %s threads for source extraction" % thread_count)
    # Assign each ND2 file to a thread, which converts it to a "fits" file
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
        # Set up a worker for each ND2 file like before
        worker_pool = Pool(thread_count)
        bfiles = [base_file for nd2_filename in image_files.filenames for base_file in base_files(nd2_filename)]
        worker_pool.map_async(source_extract, bfiles).get(timeout=sys.maxint)
        log.info("Done with Source Extractor! Took %s seconds" % round(time.time() - start, 0))
