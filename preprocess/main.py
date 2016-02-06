"""
This replaces:
    sextractor_all_directories.sh
    sextractor_directories.sh
    sextractor_txt_file.sh


It might eventually replace other files involved in making raw data available for high-level analysis.

Assumptions:
  The user has a directory with an ND2 file and some raw data from the sequencer. They will be named something like:

  15-11-18_SA15243_Cascade-TA_1nM-007.nd2
  SA15243/

  That's it!

Goal:
  Be able to run a command like "chimp preprocess" in the directory with the ND2 and NGS files and it does everything.

"""
import images
from xyz import XYZFile
from nd2reader import Nd2
import os
import subprocess
from source_extractor import SEConfig
import time
from multiprocessing import Pool
import logging
import multiprocessing

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
        log.info("fits n: %s" % n)
        subprocess.call(['fitsify', '%s%s%s.xyz' % (nd2_filename, os.path.sep, n), '%s%s%s.fits' % (nd2_filename, os.path.sep, n), '1', '2', '3'])
    log.info("Done creating fits files for %s" % nd2_filename)


def base_files(nd2_filename):
    return ["%s%s%s" % (nd2_filename, os.path.sep, os.path.splitext(filename)[0]) for filename in os.listdir(nd2_filename) if filename.endswith(".xyz")]


def source_extract(base_file):
    command = 'sextractor {base_file}.fits -PARAMETERS_NAME spot.param -CATALOG_NAME {base_file}.cat -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME {base_file}.model'
    # Don't print any output
    with open('/dev/null', 'w') as devnull:
        command = command.format(base_file=base_file).split(' ')
        subprocess.call(command)


def run():
    filenames = [nd2_filename for nd2_filename in images.get_nd2_filenames()]
    # Try to use one core per file, but top out at the number of cores that the machine has. This hasn't been proven to be optimal.
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
