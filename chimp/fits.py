import functools
import logging
from model.xyz import XYZFile
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


def get_nd2_filenames(nd2_directory):
    """
    Finds the ND2s in the current directory and returns their names without the .nd2 extension.

    """
    for f in os.listdir(nd2_directory):
        if f.endswith(".nd2"):
            yield os.path.splitext(f)[0]


def make_image_data_directory(data_directory, nd2_filename):
    """
    Creates a directory based on the ND2 filenames in order to store data derived from them.

    """
    new_directory = os.path.join(data_directory, nd2_filename)
    if not os.path.isdir(new_directory):
        os.mkdir(new_directory)


def base_files(nd2_filename):
    return ["%s%s%s" % (nd2_filename, os.path.sep, os.path.splitext(filename)[0])
            for filename in os.listdir(nd2_filename) if filename.endswith(".xyz")]


def source_extract(base_file):
    command = 'sextractor {base_file}.fits -PARAMETERS_NAME spot.param -CATALOG_NAME {base_file}.cat -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME {base_file}.model'
    # Don't print any output
    with open('/dev/null', 'w') as devnull:
        command = command.format(base_file=base_file).split(' ')
        subprocess.call(command, stdout=devnull, stderr=devnull)


def create_fits_files(data_directory, nd2_filename):
    log.info("Creating fits files for %s..." % nd2_filename)
    nd2 = Nd2(nd2_filename + ".nd2")
    indexes = []
    for n, image in enumerate(nd2):
        xyz_file = XYZFile(image)
        xyz_path = "%s.xyz" % os.path.join(data_directory, nd2_filename, n)
        with open(xyz_path, "w+") as f:
            f.write(str(xyz_file))
            indexes.append(n)
    for n in indexes:
        xyz_path = "%s.xyz" % os.path.join(data_directory, nd2_filename, n)
        subprocess.call(['fitsify', xyz_path, '%s.fits' % os.path.join(data_directory, nd2_filename, n), '1', '2', '3'])
    log.info("Done creating fits files for %s" % nd2_filename)


def run(data_directory):
    filenames = []
    for nd2_filename in get_nd2_filenames(data_directory):
        filenames.append(nd2_filename)
        make_image_data_directory(data_directory, nd2_filename)
    # Try to use one core per file, but top out at the number of cores that the machine has.
    # This hasn't been proven to be optimal.
    thread_count = min(len(filenames), multiprocessing.cpu_count())
    # Assign each ND2 file to a thread, which converts it to a "fits" file
    worker_pool = Pool(processes=thread_count)
    # The multiprocessing thing only takes an iterable with no arguments, so we use a partial function to pass
    # the directory where the files should be written
    fits_func = functools.partial(create_fits_files, data_directory)
    # KeyboardInterrupt won't behave as expected while multiprocessing unless you specify a timeout.
    # We don't want one really, so we just use the largest possible integer instead
    results = worker_pool.map_async(fits_func, filenames).get(timeout=sys.maxint)
    # Wait for the work to be finished and track how long it takes
    log.info("Starting fits file conversions.")
    start = time.time()
    results.wait()
    # log.info("results success: %s" % results.successful())
    log.info("Done with fits file conversions. Elapsed time: %s seconds" % round(time.time() - start, 0))

    # Now run source extractor (astronomy software) to do...something
    with SEConfig():
        log.info("Starting Source Extractor...")
        start = time.time()
        # Set up a worker for each ND2 file like before
        worker_pool = Pool(thread_count)
        files = [base_file for nd2_filename in filenames for base_file in base_files(nd2_filename)]
        results = worker_pool.map_async(source_extract, files).get(timeout=sys.maxint)

        # Wait for the work to be done
        results.wait()
    log.info("Done with Source Extractor! Took %s seconds" % round(time.time() - start, 0))
