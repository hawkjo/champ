from copy import deepcopy
import numpy as np
import misc


class FastqTileRCs(object):
    """ A class for fastq tile coordinates. """
    def __init__(self, key, read_names):
        self.key = key
        self.rcs = np.array([map(int, name.split(':')[-2:]) for name in read_names])
        self._image = None

    def set_fastq_image_data(self, offset, scale, scaled_dims, width):
        self.offset = offset
        self.scale = scale
        self.image_shape = scaled_dims
        self.width = width  # width in um
        self.mapped_rcs = scale * (self.rcs + np.tile(offset, (self.rcs.shape[0], 1)))

    def rotate_data(self, degrees):
        self.mapped_rcs = np.dot(self.mapped_rcs, misc.right_rotation_matrix(degrees))
        self.mapped_rcs -= np.tile(self.mapped_rcs.min(axis=0), (self.mapped_rcs.shape[0], 1))
        self.image_shape = self.mapped_rcs.max(axis=0) + 1
        return self.image_shape

    @property
    def image(self):
        if self._image is None:
            self._image = np.zeros(self.image_shape)
            self._image[self.mapped_rcs.astype(np.int)[:, 0], self.mapped_rcs.astype(np.int)[:, 1]] = 1
        return self._image

    def fft_align_with_im(self, image_data):
        im_data_im_shapes = set(a.shape for a in image_data.all_ffts.values())
        assert len(im_data_im_shapes) <= 2, im_data_im_shapes

        # Make the ffts
        fq_im_fft_given_shape = {}
        for shape in im_data_im_shapes:
            padded_fq_im = misc.pad_to_size(self.image, shape)
            fq_im_fft_given_shape[shape] = np.fft.fft2(padded_fq_im)

        # Align
        self.best_max_corr = float('-inf')
        for im_key, im_data_fft in image_data.all_ffts.items():
            fq_im_fft = fq_im_fft_given_shape[im_data_fft.shape]
            cross_corr = abs(np.fft.ifft2(np.conj(fq_im_fft) * im_data_fft))
            max_corr = cross_corr.max()
            max_idx = misc.max_2d_idx(cross_corr)

            if max_corr > self.best_max_corr:
                self.best_im_key = im_key
                self.best_max_corr = max_corr
                self.align_tr = np.array(max_idx) - self.image.shape
        return self.key, self.best_im_key, self.best_max_corr, self.align_tr

    def get_new_aligned_rcs(self, new_fq_w=None, new_degree_rot=0, new_tr=(0, 0)):
        """Returns aligned rcs. Only works when image need not be flipped or rotated."""
        if new_fq_w is None:
            new_fq_w = self.width
        aligned_rcs = deepcopy(self.mapped_rcs)
        aligned_rcs = np.dot(aligned_rcs, misc.right_rotation_matrix(new_degree_rot))
        aligned_rcs -= np.tile(aligned_rcs.min(axis=0), (aligned_rcs.shape[0], 1))
        aligned_rcs *= float(new_fq_w) / self.width
        aligned_rcs += np.tile(self.align_tr + new_tr, (aligned_rcs.shape[0], 1))
        return aligned_rcs

    def set_aligned_rcs(self):
        self.aligned_rcs = self.get_new_aligned_rcs()

    def set_aligned_rcs_given_transform(self, lbda, theta, offset):
        """Performs transform calculated in FastqImageAligner.least_squares_mapping."""
        A = np.zeros((2*len(self.rcs), 4))
        for i, pt in enumerate(self.rcs):
            xir, yir = pt
            A[2*i, :] = [xir, -yir, 1, 0]
            A[2*i+1, :] = [yir,  xir, 0, 1]

        x = np.array([lbda * np.cos(theta),
                      lbda * np.sin(theta),
                      offset[0],
                      offset[1]])

        # First update w since it depends on previous scale setting
        self.width = lbda * float(self.width) / self.scale
        self.scale = lbda
        self.rotation = theta
        self.offset = offset
        self.aligned_rcs = np.dot(A, x).reshape((len(self.rcs), 2))

    def set_correlation(self, im):
        """Sets alignment correlation. Only works when image need not be flipped or rotated."""
        self.best_max_corr = sum(im[pt[0], pt[1]] for pt in self.aligned_rcs
                                 if 0 <= pt[0] < im.shape[0] and 0 <= pt[1] < im.shape[1])

    def set_snr(self, snr):
        self.snr = snr

    def set_snr_with_control_corr(self, control_corr):
        self.snr = self.best_max_corr / control_corr
