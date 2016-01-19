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
from pathos.multiprocessing import ProcessPool
import logging
import multiprocessing

log = logging.getLogger()
log.addHandler(logging.StreamHandler())
log.setLevel(logging.DEBUG)


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
        subprocess.call(['fitsify', '%s%s%s.xyz' % (nd2_filename, os.path.sep, n), '%s%s%s.fits' % (nd2_filename, os.path.sep, n), '1', '2', '3'])
    log.info("Done creating fits files for %s" % nd2_filename)


if __name__ == "__main__":
    filenames = [nd2_filename for nd2_filename in images.get_nd2_filenames()]
    filecount = len(filenames)
    cores = multiprocessing.cpu_count()
    p = ProcessPool(nodes=min(filecount, cores))
    results = p.amap(create_fits_files, filenames)
    while not results.ready():
        time.sleep(3)
        log.debug(time.time())
    log.info("Done with all fits file conversions")
    with SEConfig():
        log.info("Starting Source Extractor...")
    log.info("Done with Source Extractor!")

"""
Now do:

# Run sextractor
sextractor $base.fits -PARAMETERS_NAME spot.param -CATALOG_NAME $base.cat -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME $base.model

echo "****** clean up"
# Clean up
if [[ $CLEAN -gt 0 ]] ; then
  rm $base.xyz default.conv spot.param default.sex
fi

if [[ $CLEAN -gt 1 ]] ; then
  rm $base.model
fi

if [[ $CLEAN -gt 2 ]] ; then
  rm $base.fits
fi

exit 0
"""