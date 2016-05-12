from chimp import convert
import logging
import os

log = logging.getLogger(__name__)


def main(clargs):
    print(clargs.flipud, clargs.fliplr)
    exit()
    for directory in clargs.tif_directories:
        if not os.path.isdir(directory):
            log.debug("Skipping non-directory %s" % directory)
            continue
        tif_filenames = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".ome.tif")]
        if tif_filenames:
            hdf5_file_path = "%s.h5" % directory
            log.debug("About to convert %d files in %s to HDF5 file: %s" % (len(tif_filenames),
                                                                            directory,
                                                                            hdf5_file_path))
            convert.tif_dir_to_hdf5(hdf5_file_path, tif_filenames, clargs.flipud, clargs.fliplr)
