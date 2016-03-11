import os
import numpy as np
from misc import median_normalize as normalize_median


class ImageData(object):
    """ A class for image data to be correlated with fastq coordinate data. """
    def __init__(self, fpath=None, objective=None, im=None, fname='<name unknown>', median_normalize=False):
        if im is not None:
            if fpath:
                self.fpath = fpath
                fname = fpath.split('/')[-1]
            self.set_im_from_ndarray(im, objective, fname, median_normalize)
        elif fpath:
            self.set_im_from_file(fpath, objective, median_normalize)

    def set_im_from_ndarray(self, im, objective, fname='<name unknown>', median_normalize=False):
        assert isinstance(im, np.ndarray), 'Image not numpy ndarray'
        self.im = im
        self.set_objective(objective)
        self.fname = fname
        self.bname = os.path.splitext(self.fname)[0]
        if median_normalize:
            normalize_median(self.im)

    def set_im_from_file(self, fpath, objective, median_normalize=False):
        self.fpath = fpath
        self.fname = fpath.split('/')[-1]
        self.bname, ext = os.path.splitext(self.fname)
        if ext == '.npy':
            self.im = np.load(self.fpath)
        else:
            raise ValueError('Image type not accepted: %s' % self.fname)
        self.set_objective(objective)
        if median_normalize:
            normalize_median(self.im)

    def set_objective(self, objective):
        assert objective in [20, 40, 60], 'Accepted objectives are 20, 40, and 60'
        self.objective = objective
        self.um_per_pixel = 16.0 / self.objective
        self.um_dims = self.um_per_pixel * np.array(self.im.shape)
