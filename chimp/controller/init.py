from chimp import convert
import logging

log = logging.getLogger(__name__)


def main(clargs):
    paths = convert.get_all_tif_paths(clargs.image_directory)
    # directories will have ".h5" appended to them to come up with the HDF5 names
    # tifs are relative paths to each tif file
    convert.main(paths, clargs.flipud, clargs.fliplr)
