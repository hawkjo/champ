import numpy as np
from pathos.multiprocessing import ProcessingPool
from misc import next_power_of_2


class ImageData(object):
    """A class for image data to be correlated with fastq coordinate data."""
    def __init__(self, im, objective, im_type='np'):
        if im_type in ['np', 'numpy']:
            assert isinstance(im, np.ndarray), 'Image not numpy ndarray'
            self.im = im
        elif im_type == 'tif':
            self.im = np.array(Image.open(im))
        else:
            raise ValueError('Image type not accepted: %s' % im_type)

        assert objective in set([40, 60]), 'Accepted objectives are 40 and 60'
        self.objective = objective

        self.um_per_pixel = 16.0 / self.objective
        self.um_dims = self.um_per_pixel * np.array(self.im.shape)

    def D4_im_given_idx(self, idx):
        flip, rot = idx
        if flip:
            flip_im = np.fliplr(self.im)
        else:
            flip_im = self.im
        return np.rot90(flip_im, k=(rot/90))

    def iterate_D4_images(self):
        for idx in self.iterate_D4_idxs():
            yield idx, self.D4_im_given_idx(idx)

    def iterate_D4_idxs(self):
        for flip in [0, 1]:
            for rot in [0, 90, 180, 270]:
                yield (flip, rot)

    def D4_ffts(self, padding=(0, 0), processors=1, force=False):
        """Makes images and ffts of all flips and 90 degree rotations (i.e. D4)"""
        if hasattr(self, 'all_ffts') and self.all_ffts and not force:
            return
        self.fft_padding = padding
        pool = ProcessingPool(processors)
        ffts = pool.map(self.single_fft, self.iterate_D4_idxs())
        self.all_ffts = {idx: fft for idx, fft in zip(self.iterate_D4_idxs(), ffts)}
        del pool

    def set_single_fft(self, idx, padding):
        self.fft_padding = padding
        self.all_ffts = {idx: self.single_fft(idx)}

    def single_fft(self, idx):
        im = self.D4_im_given_idx(idx)
        totalx, totaly = np.array(self.fft_padding)+np.array(im.shape)
        w = next_power_of_2(totalx)
        h = next_power_of_2(totaly)
        padded_im = np.pad(im,
                           ((self.fft_padding[0], w-totalx), (self.fft_padding[1], h-totaly)),
                           mode='constant')
        return np.fft.fft2(padded_im)
