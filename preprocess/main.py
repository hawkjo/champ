import images
import logging
import multiprocessing
from multiprocessing import Pool
from nd2reader import Nd2
import os
from source_extractor import SEConfig
import subprocess
import time
from xyz import XYZFile


log = logging.getLogger(__name__)


def create_fits_files(nd2_filename):
    log.info("Creating fits files for %s..." % nd2_filename)
    nd2 = Nd2(nd2_filename + ".nd2")
    images.make_image_data_directory(nd2_filename)
    indexes = []
    for n, image in enumerate(nd2):
        xyz_file = XYZFile(image)
        with open("%s%s%s.xyz" % (nd2_filename, os.path.sep, n), "w+") as f:
            f.write(str(xyz_file))
            indexes.append(n)
    for n in indexes:
        subprocess.call(['fitsify', '%s%s%s.xyz' % (nd2_filename, os.path.sep, n), '%s%s%s.fits' % (nd2_filename, os.path.sep, n), '1', '2', '3'])
    log.info("Done creating fits files for %s" % nd2_filename)


def base_files(nd2_filename):
    return ["%s%s%s" % (nd2_filename, os.path.sep, os.path.splitext(filename)[0]) for filename in os.listdir(nd2_filename) if filename.endswith(".xyz")]


def source_extract(base_file):
    command = 'sextractor {base_file}.fits -PARAMETERS_NAME spot.param -CATALOG_NAME {base_file}.cat -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME {base_file}.model'
    # Don't print any output
    with open('/dev/null', 'w') as devnull:
        command = command.format(base_file=base_file).split(' ')
        subprocess.call(command, stdout=devnull, stderr=devnull)


def run():
    filenames = [nd2_filename for nd2_filename in images.get_nd2_filenames()]
    # Try to use one core per file, but top out at the number of cores that the machine has.
    # This hasn't been proven to be optimal.
    thread_count = min(len(filenames), multiprocessing.cpu_count())

    # Assign each ND2 file to a thread, which converts it to a "fits" file
    worker_pool = Pool(processes=thread_count)
    results = worker_pool.map_async(create_fits_files, filenames)

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
        results = worker_pool.map_async(source_extract, [base_file for nd2_filename in filenames for base_file in base_files(nd2_filename)])

        # Wait for the work to be done
        results.wait()
    log.info("Done with Source Extractor! Took %s seconds" % round(time.time() - start, 0))
