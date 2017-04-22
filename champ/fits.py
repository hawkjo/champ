from astropy.io import fits
from champ.grid import GridImages
import h5py
import logging
import numpy as np
import os
from scipy import ndimage
from skimage.filters import threshold_otsu
import subprocess
import threading
from Queue import Queue, Empty
import time
import sqlite3


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
            for filename in os.listdir(h5_filename) if filename.endswith(".fits")]


def ensure_image_data_directory_exists(h5_filename):
    """
    Creates a directory based on the HDF5 filenames in order to store data derived from them.

    """
    new_directory = os.path.join(h5_filename)
    if not os.path.isdir(new_directory):
        os.mkdir(new_directory)


def otsu_cluster_func(h5_filename, condition, results_queue):
    log.info("Finding Otsu clusters for %s" % condition)
    with h5py.File(h5_filename, 'r') as h5:
        images = h5[condition]
        for channel in images.keys():
            grid = GridImages(images, channel)
            for image in grid:
                center_of_masses = find_center_of_masses_using_otsu(image)
                results_queue.put((center_of_masses, condition, image.row, image.column))


def find_center_of_masses_using_otsu(image):
    threshold = threshold_otsu(image)
    mask_pixels = (image > threshold)
    mask = ndimage.binary_closing(ndimage.binary_opening(mask_pixels))
    label_image, num_labels = ndimage.label(mask)
    center_of_masses = ndimage.center_of_mass(image, label_image, range(num_labels + 1))
    return center_of_masses


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


def source_extract(base_file):
    command = '/usr/bin/sextractor {base_file}.fits -PARAMETERS_NAME spot.param -CATALOG_NAME STDOUT -CATALOG_TYPE ASCII -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME {base_file}.model -VERBOSE_TYPE QUIET'
    # Don't print any output
    command = command.format(base_file=base_file).split(' ')
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cluster_data, errors = process.communicate()
    if errors:
        log.error("Error using Source Extractor on file: %s" % base_file)
    with open("{base_file}.clusters.se".format(base_file=base_file), 'w') as f:
        f.write(cluster_data)


def create_fits_files(h5_filename, condition):
    log.info("Creating fits files for %s" % h5_filename)
    h5 = h5py.File(h5_filename)
    for channel in h5[condition].keys():
        grid = GridImages(h5[condition], channel)
        for image in grid:
            fits_path = '%s.fits' % os.path.join('%s.%s' % (condition, image.index))
            # Source Extractor can handle at most 32-bit values, so we have to cast down from our 64-bit images or
            # else it will throw an error. We clip to ensure there's no overflow, although this seems improbable
            # given that most cameras are 16 bit
            clipped_image = np.clip(image, 0, 2**32-1).astype(np.uint32)
            hdu = fits.PrimaryHDU(clipped_image)
            hdu.writeto(fits_path, clobber=True)
    log.info("Done creating fits files for %s" % condition)


def main(image_directory):
    db = sqlite3.connect(os.path.join(image_directory, 'champ.db'))
    cursor = db.cursor()
    cursor.execute("DROP TABLE IF EXISTS clusters")
    cursor.execute("DROP TABLE IF EXISTS fields_of_view")
    cursor.execute("CREATE TABLE clusters ("
                   "field_of_view_id INTEGER, "
                   "r FLOAT, "
                   "c FLOAT, PRIMARY KEY (field_of_view_id, r, c)) WITHOUT ROWID")
    cursor.execute("CREATE TABLE fields_of_view ("
                   "id INTEGER PRIMARY KEY ASC, "
                   "condition MEDIUMINT UNSIGNED, "
                   "r TINYINT UNSIGNED, "
                   "c SMALLINT UNSIGNED)")
    db.commit()

    h5_filename = os.path.join(image_directory, 'images.h5')
    with h5py.File(h5_filename, 'r') as h5:
        conditions = h5.keys()
    # Assign each HDF5 file to a thread, which converts it to a "fits" file
    results_queue = Queue()
    threads = []
    for condition in conditions:
        thread = threading.Thread(target=thread_find_clusters_otsu, args=(h5_filename, condition, results_queue))
        thread.start()
        threads.append(thread)
    for thread in threads:
        thread.join()
    cursor = db.cursor()
    while True:
        try:
            center_of_masses, condition, fov_row, fov_column = results_queue.get_nowait()
            cursor.execute("INSERT INTO fields_of_view (condition, r, c) VALUES ({condition}, {row}, {column})".format(
                condition=condition,
                row=fov_row,
                column=fov_column
            ))
            # write_field_of_view(condition, fov_row, fov_column)
            # write_cluster_locations(center_of_masses, cursor.lastrowid, cursor)
            query = "INSERT INTO clusters VALUES ({field_of_view_id}, ?, ?)".format(field_of_view_id=cursor.lastrowid)
            cursor.executemany(query, center_of_masses)
        except Empty:
            break
    db.commit()
    db.close()


def thread_find_clusters_otsu(h5_filename, condition, results_queue):
    # Find clusters with Otsu thresholding
    start = time.time()
    otsu_cluster_func(h5_filename, condition, results_queue)
    # otsu_partial_func = functools.partial(otsu_cluster_func, image_file)
    # worker_pool.map_async(otsu_partial_func, conditions).get(timeout=sys.maxint)
    log.info("Done with cluster location. Elapsed time: %s seconds" % round(time.time() - start, 0))


# def find_clusters_source_extractor(worker_pool, image_file, conditions):
#     # Find clusters with Source Extractor
#     log.info("Starting fits file conversions.")
#     # The multiprocessing thing only takes an iterable with no arguments, so we use a partial function to pass
#     # the directory where the files should be written
#     fits_func = functools.partial(create_fits_files, image_file)
#     # KeyboardInterrupt won't behave as expected while multiprocessing unless you specify a timeout.
#     # We don't want one really, so we just use the largest possible integer instead
#     start = time.time()
#     worker_pool.map_async(fits_func, image_files.directories).get(timeout=sys.maxint)
#     log.info("Done with fits file conversions. Elapsed time: %s seconds" % round(time.time() - start, 0))
#
#     # Now run source extractor to find the coordinates of points
#     with SEConfig():
#         log.info("Starting Source Extractor...")
#         start = time.time()
#         # Set up a worker for each HDF5 file like before
#         base_files = [base_file for h5_filename in image_files.directories
#                       for base_file in get_base_file_names(h5_filename)]
#         worker_pool.map_async(source_extract, base_files).get(timeout=sys.maxint)
#         log.info("Done with Source Extractor! Took %s seconds" % round(time.time() - start, 0))
