import logging
from champ import convert, initialize

log = logging.getLogger(__name__)


def main(clargs):
    metadata = initialize.load_metadata(clargs.image_directory)
    # directories will have ".h5" appended to them to come up with the HDF5 names
    # tifs are relative paths to each tif file
    log.debug("Preprocessing images.")
    paths = convert.get_all_tif_paths(clargs.image_directory)
    log.debug("About to convert TIFs to HDF5.")
    convert.main(paths, metadata['flipud'], metadata['fliplr'], clargs.min_column, clargs.max_column)
    log.debug("Done converting TIFs to HDF5.")
