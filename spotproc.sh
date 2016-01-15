#!/bin/bash
# Process a spot text file through source extractor.
#
# Assumes fitsify ("make fitsify") 
#         sextractor ("apt-get install sextractor")
#
# Syntax: spotproc basename
#   E.g.  spotproc 150403_M01012_210
#
# will process a text 512x512 file 150403_M01012_210.txt 
#
# Output includes:
#
#   150403_M01012_210.fits  = FITS format image
#   150403_M01012_210.cat   = Catalog output from sextractor
#   150403_M01012_210.model = Sextractor model of image for checking
#

base=$1
shift 1

CLEAN=3

eval $@

# Convert number<tab>number<tab>number... to x,y coords and number
tr '\t' '\n' < $base.txt | awk '{printf "%3d %3d %5d\n", (NR-1)%512, int((NR-1)/512), $1}' > $base.xyz
#tr '\t' '\n' < $base.txt | awk '{printf "%3d %3d %e\n", (NR-1)%512, int((NR-1)/512), $1}' > $base.xyz

# Convert to a fits file
fitsify $base.xyz $base.fits 1 2 3

echo "****** create config files"
# Create a sextractor default configuration file?  Nah, use specialized.
# sextractor -d > default.sex
# Changed my mind.
cat > default.sex <<EOF
DETECT_THRESH 2
DEBLEND_NTHRESH 64
DEBLEND_MINCONT 0.00005
EOF


# Create a parameter file to tell sextractor what output we want
cat > spot.param <<EOF
X_IMAGE
Y_IMAGE
FLUX_AUTO
FLUXERR_AUTO
FLAGS
A_IMAGE
B_IMAGE
THETA_IMAGE
EOF

# Create a default convolution file for extractor
cat > default.conv <<EOF
CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
EOF

echo "****** running sextractor"
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
