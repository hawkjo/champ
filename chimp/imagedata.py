import numpy as np
from chimp import misc


class ImageData(object):
    """A class for image data to be correlated with fastq coordinate data."""
    def __init__(self, fname, objective, image):
        assert isinstance(image, np.ndarray), 'Image not numpy ndarray'
        assert objective in [20, 40, 60], 'Accepted objectives are 20, 40, and 60'
        self.fname = str(fname)
        self.fft = None
        self.image = image
        self.median_normalize()
        self.objective = objective
        self.um_per_pixel = 16.0 / self.objective
        self.um_dims = self.um_per_pixel * np.array(self.image.shape)

    def median_normalize(self):
        med = np.median(self.image)
        self.image = self.image.astype('float', copy=False, casting='safe')
        self.image /= float(med)
        self.image -= 1.0

    def set_fft(self, padding):
        totalx, totaly = np.array(padding) + np.array(self.image.shape)
        w = misc.next_power_of_2(totalx)
        h = misc.next_power_of_2(totaly)
        padded_im = np.pad(self.image,
                           ((int(padding[0]), int(w-totalx)), (int(padding[1]), int(h-totaly))),
                           mode='constant')
        self.fft = np.fft.fft2(padded_im)
