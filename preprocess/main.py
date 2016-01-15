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
from preprocess import images
from preprocess.xyz import XYZFile
from nd2reader import Nd2
import os
import subprocess
from source_extractor import SEConfig

if __name__ == "__main__":
    for nd2_filename in images.get_nd2_filenames():
        with SEConfig(nd2_filename) as sec:
            nd2 = Nd2(nd2_filename + ".nd2")
            images.make_image_data_directory(nd2_filename)
            os.chdir(nd2_filename)
            for n, image in enumerate(nd2):
                xyz_file = XYZFile(image)
                with open("%s.xyz" % n, "w+") as f:
                    f.write(str(xyz_file))
                subprocess.call(['fitsify', '%s.xyz' % n, '%s.fits' % n, '1', '2', '3'])
            os.chdir('..')

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