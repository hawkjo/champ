import h5py
import numpy as np
import os
import sys
import hdf5_tools
from skimage.filters import threshold_otsu
from scipy import ndimage


def thresh_and_locate_spots(h5_fpath):
    out_dir = os.path.splitext(h5_fpath)[0] + '_sex_files'
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    with h5py.File(h5_fpath, 'r') as f:
        for channel, group in f.items():
            for dset_name, dset in group.items():
                sys.stdout.write('.')
                sys.stdout.flush()

                out_fname = hdf5_tools.bname_given_channel_and_dset_name(channel, dset_name) + '.cat'
                out_fpath = os.path.join(out_dir, out_fname)

                im = np.array(dset)
                thresh = threshold_otsu(im)
                mask = (im > thresh)
                mask = ndimage.binary_closing(ndimage.binary_opening(mask))
                label_im, num_labels = ndimage.label(mask)
                c_of_masses = ndimage.center_of_mass(im, label_im, range(num_labels+1))

                write_locs_in_sextractor_format(c_of_masses, out_fpath)


def write_locs_in_sextractor_format(locs, out_fpath):
    with open(out_fpath, 'w') as out:
        out.write("""#   1 X_IMAGE                Object position along x                                    [pixel]
#   2 Y_IMAGE                Object position along y                                    [pixel]
#   3 FLUX_AUTO              Flux within a Kron-like elliptical aperture                [count]
#   4 FLUXERR_AUTO           RMS error for AUTO flux                                    [count]
#   5 FLAGS                  Extraction flags                                          
#   6 A_IMAGE                Profile RMS along major axis                               [pixel]
#   7 B_IMAGE                Profile RMS along minor axis                               [pixel]
#   8 THETA_IMAGE            Position angle (CCW/x)                                     [deg]
""")
        out.write(
            '\n'.join(
                '\t'.join(map(str, [c+1, r+1] + [0]*6))  # Sextrator does 1-based x and y
                for r, c in locs
            )
        )


if __name__ == '__main__':
    usage_fmt = '{} <h5_fpath>'.format(sys.argv[0])
    if len(sys.argv) != len(usage_fmt.split()):
        sys.exit(usage_fmt)
    sextractor_hdf5(sys.argv[1])
