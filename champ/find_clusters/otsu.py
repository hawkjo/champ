from champ.grid import GridImages
import logging
import os
from threading import Thread
from Queue import Queue
import h5py
from skimage.filters import threshold_otsu
from scipy import ndimage


log = logging.getLogger(__name__)


def main(image_files, thread_count):
    queue = Queue()
    for t in range(thread_count):
        thread = Thread(target=find_clusters, args=(queue,))
        thread.daemon = True
        thread.start()
    for h5_base_name in image_files.directories:
        queue.put(h5_base_name)
    queue.join()


def find_clusters(queue):
    while True:
        h5_base_name = queue.get()
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
        queue.task_done()


def write_cluster_locations(locations, out_filepath):
    with open(out_filepath, 'w') as out:
        out.write('\n'.join("%s\t%s" % (r, c) for c, r in locations))
