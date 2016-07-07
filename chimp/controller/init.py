from chimp import convert, fits
import logging
import yaml
import os


log = logging.getLogger(__name__)


def main(clargs):
    paths = convert.get_all_tif_paths(clargs.image_directory)
    # directories will have ".h5" appended to them to come up with the HDF5 names
    # tifs are relative paths to each tif file
    # TODO: Make it save the metadata to disk - that way we don't have to pass the image directory or mapped_reads dir again
    log.debug("About to convert TIFs to HDF5.")
    convert.main(paths, clargs.flipud, clargs.fliplr)
    log.debug("Done converting TIFs to HDF5.")
    log.debug("Fitsifying images from HDF5 files.")
    fits.main(clargs.image_directory)


def save_metadata(image_directory, chip_name, mapped_reads_dir, alignment_channel,
                  microns_per_pixel, chip_type, ports_on_right):
    with open(os.path.join(image_directory, 'champ.yaml'), 'w+') as f:
        data = {'chip_name': chip_name,
                'mapped_reads': os.path.abspath(mapped_reads_dir),
                'microns_per_pixel': microns_per_pixel,
                'chip_type': chip_type,
                'ports_on_right': ports_on_right,
                'alignment_channel': alignment_channel}
        yaml.dump(data, f)
