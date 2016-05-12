from chimp.convert import tif_dir_to_hdf5
import os


def main(clargs):
    for directory in clargs.tif_directories:
        tif_filenames = [f for f in os.listdir(directory) if f.endswith(".ome.tif")]
        print(tif_filenames)
        if tif_filenames:
            hdf5_file_path = "%s.h5" % directory
            tif_dir_to_hdf5(hdf5_file_path, tif_filenames, clargs.flipud, clargs.fliplr)
