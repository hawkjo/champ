from chimp.convert import tif_dir_to_hdf5


def main(clargs):
    tif_dir_to_hdf5(clargs.hdf5_file_path, clargs.tif_file_paths, clargs.flipud, clargs.fliplr)
