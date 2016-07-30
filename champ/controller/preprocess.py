from champ import convert, fits, initialize
import logging


log = logging.getLogger(__name__)


def main(clargs):
    experiment = initialize.load(clargs.image_directory)
    paths = convert.get_all_tif_paths(clargs.image_directory)
    # directories will have ".h5" appended to them to come up with the HDF5 names
    # tifs are relative paths to each tif file
    log.debug("About to convert TIFs to HDF5.")
    convert.main(paths, experiment['flipud'], experiment['fliplr'], clargs.enhance_contrast)
    log.debug("Done converting TIFs to HDF5.")
    log.debug("Fitsifying images from HDF5 files.")
    fits.main(clargs.image_directory)
