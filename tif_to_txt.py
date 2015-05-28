import sys
import os
import glob
import numpy as np
import tifffile


def tif_to_txt(im_dir, fmt):
    fpaths = glob.glob(os.path.join(im_dir, '*.tif'))
    for fpath in fpaths:
        tif = tifffile.imread(fpath)
        if len(tif.shape) == 2:
            out_fname = fpath[:-3] + 'txt'
            np.savetxt(out_fname, tif, fmt=fmt, delimiter='\t')
        if len(tif.shape) == 3:
            bname = fpath[:-4]
            fnum_fmt = 'f%0' + str(len(str(tif.shape[0]))) + 'd'
            for i in range(tif.shape[0]):
                out_fname = bname + fnum_fmt % i + '.txt'
                np.savetxt(out_fname, tif[i], fmt=fmt, delimiter='\t')
    

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit('Usage: %s <int or float> <directory>' % sys.argv[0])

    if sys.argv[1].lower() == 'int':
        fmt = '%d'
    elif sys.argv[1].lower() == 'float':
        fmt = '%f'
    else:
        sys.exit('Format must be int or float')

    im_dir = sys.argv[2]

    tif_to_txt(im_dir, fmt)
