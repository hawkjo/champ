import numpy as np
from chimp import misc
from skimage import io


class ImageData(object):
    """A class for image data to be correlated with fastq coordinate data."""
    def __init__(self, filename, um_per_pixel, image):
        assert isinstance(image, np.ndarray), 'Image not numpy ndarray'
        self.fname = str(filename)
        self.fft = None
        self.image = image
        self.median_normalize()
        self.um_per_pixel = um_per_pixel
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
        io.imshow(padded_im)
        io.show()
        self.fft = np.fft.fft2(padded_im)
