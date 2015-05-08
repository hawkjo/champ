import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import local_config
import sextraction
from scipy.spatial import KDTree
from Bio import SeqIO
from imagedata import ImageData
from fastqtilercs import FastqTileRCs
from misc import pad_to_size, max_2d_idx


class FastqImageCorrelator(object):
    """A class to find the alignment of fastq data and image data."""
    def __init__(self, project_name):
        self.project_name = project_name
        self.fastq_tiles = {}
        self.fastq_tiles_list = []
        self.fastq_tiles_keys = []
        self.image_data = ImageData()
        self.w_fq_tile_min = 895  # um
        self.w_fq_tile_max = 937  # um
        self.w_fq_tile = 937  # um

    def load_phiX(self):
        for key, tile in local_config.phiX_rcs_given_project_name(self.project_name).items():
            self.fastq_tiles[key] = FastqTileRCs(key, tile)
        self.fastq_tiles_keys = [key for key, tile in sorted(self.fastq_tiles.items())]
        self.fastq_tiles_list = [tile for key, tile in sorted(self.fastq_tiles.items())]

    def set_image_data(self, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], ImageData):
            self.image_data = args[0]
        else:
            self.image_data = ImageData(*args, **kwargs)

    def set_fastq_tile_mappings(self):
        """Calculate parameters for mapping fastq tiles for ffts."""
        assert self.image_data is not None, 'No image data loaded.'
        assert self.fastq_tiles != {}, 'No fastq data loaded.'

        self.all_data = np.concatenate([tile.rcs for key, tile in self.fastq_tiles.items()])
    
        x_min, y_min = self.all_data.min(axis=0)
        x_max, y_max = self.all_data.max(axis=0)
    
        self.fq_im_offset = np.array([-x_min, -y_min])
        self.fq_im_scale = (float(self.w_fq_tile) / (x_max-x_min)) / self.image_data.um_per_pixel
        self.fq_im_scaled_maxes = self.fq_im_scale * np.array([x_max-x_min, y_max-y_min])
        self.fq_im_scaled_dims = (self.fq_im_scaled_maxes + [1, 1]).astype(np.int)

    def set_all_fastq_image_data(self, verbose=True):
        for key, tile in self.fastq_tiles.items():
            tile.set_fastq_image_data(self.fq_im_offset,
                                      self.fq_im_scale,
                                      self.fq_im_scaled_dims,
                                      self.w_fq_tile,
                                      verbose=verbose)

    def rotate_all_fastq_data(self, degrees):
        im_shapes = []
        for tile in self.fastq_tiles_list:
            im_shapes.append(tile.rotate_data(degrees))
        self.fq_im_scaled_dims = np.array(im_shapes).max(axis=0)
        for tile in self.fastq_tiles_list:
            tile.image_shape = self.fq_im_scaled_dims


    def imreg_align(self):
        for key, tile in sorted(self.fastq_tiles.items()):
            tile.imreg_align_with_im(self.image_data.im)

    def fft_align_tile(self, tile):
        return tile.fft_align_with_im(self.image_data)

    def fft_align_tile_with_im(self, tile):
        im_data_im_shapes = set(a.shape for a in self.image_data.all_ffts.values())
        assert len(im_data_im_shapes) <= 2, im_data_im_shapes

        # Make the ffts
        fq_image = tile.image()
        fq_im_fft_given_shape = {}
        for shape in im_data_im_shapes:
            padded_fq_im = pad_to_size(fq_image, shape)
            fq_im_fft_given_shape[shape] = np.fft.fft2(padded_fq_im)

        # Align
        best_max_corr = float('-inf')
        for im_key, im_data_fft in self.image_data.all_ffts.items():
            fq_im_fft = fq_im_fft_given_shape[im_data_fft.shape]
            cross_corr = abs(np.fft.ifft2(np.conj(fq_im_fft) * im_data_fft))
            max_corr = cross_corr.max()
            max_idx = max_2d_idx(cross_corr)

            if max_corr > best_max_corr:
                best_im_key = im_key
                best_max_corr = max_corr
                align_tr = np.array(max_idx) - fq_image.shape
        #print 'Result:', tile.key, best_im_key, best_max_corr, align_tr
        return tile.key, best_im_key, best_max_corr, align_tr

    def fft_align(self, processors, recalc_fft=True, verbose=True):
        if verbose:
            print 'Set fastq tile mappings'
        self.set_fastq_tile_mappings()
        if verbose:
            print 'Image D4 ffts'
        self.image_data.D4_ffts(padding=self.fq_im_scaled_dims,
                                processors=processors,
                                force=recalc_fft)
        if verbose:
            print 'Fastq images and ffts'
        self.set_all_fastq_image_data(verbose=True)
        if verbose:
            print 'Aligning'
        self.best_corr = 0
        self.max_corrs = []
        for tile in self.fastq_tiles_list:
            fq_key, im_key, max_corr, align_tr = self.fft_align_tile(tile)
            self.max_corrs.append(max_corr)
            if max_corr > self.best_corr:
                self.best_corr = max_corr
                self.best_fq_key = fq_key
                self.best_im_key = im_key
                self.best_align_tr = align_tr
        if verbose:
            print 'Best result:', self.best_corr, self.best_fq_key, self.best_im_key, self.best_align_tr

    def show_alignment(self, fq_key, im_key, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        im = self.image_data.D4_im_given_idx(im_key)
        ax.imshow(im, cmap=plt.get_cmap('Blues'), label='Ilya')
        v = plt.axis()
        fq_tile = self.fastq_tiles[fq_key]
        fq_im = fq_tile.image()[-fq_tile.align_tr[0]:, -fq_tile.align_tr[1]:]
        fq_im = ndimage.filters.gaussian_filter(fq_im, 3)
        fq_im[fq_im < 0.01] = None
        plt.imshow(fq_im, cmap=plt.get_cmap('Reds'), alpha=0.5, label='Fastq')
        #aligned_rcs = self.fastq_tiles[fq_key].aligned_rcs()
        #ax.plot(aligned_rcs[:, 0], aligned_rcs[:, 1], 'r.', alpha=0.5)
        plt.axis(v)

    def set_sexcat(self, sexcat):
        assert isinstance(sexcat, sextraction.Sextraction)
        self.sexcat = sexcat

    def set_sexcat_from_file(self, fpath):
        self.sexcat = sextraction.Sextraction(fpath)

    def find_hits(self, tile_key):
        # First, restrict to points in the frame, keeping rcs along for backtracking
        self.aligned_rcs_in_frame = []
        self.rcs_in_frame = []
        rcs = self.fastq_tiles[tile_key].rcs.astype(np.int)
        im_shape = self.image_data.im.shape
        for i, pt in enumerate(self.fastq_tiles[tile_key].aligned_rcs):
            if 0 <= pt[0] < im_shape[0] and 0 <= pt[1] < im_shape[1]:
                self.aligned_rcs_in_frame.append(pt)
                self.rcs_in_frame.append(rcs[i])

        # Next find nearest neighbors
        sexcat_tree = KDTree(self.sexcat.point_rcs)
        aligned_tree = KDTree(self.aligned_rcs_in_frame)

        # All indices are in the order (sexcat_idx, aligned_in_frame_idx)
        sexcat_to_aligned_idxs = set()
        for i, pt in enumerate(self.sexcat.point_rcs):
            dist, idx = aligned_tree.query(pt)
            sexcat_to_aligned_idxs.add((i, idx))

        aligned_to_sexcat_idxs_rev = set()
        for i, pt in enumerate(self.aligned_rcs_in_frame):
            dist, idx = sexcat_tree.query(pt)
            aligned_to_sexcat_idxs_rev.add((idx, i))

        # Find categories of hits
        mutual_hits = sexcat_to_aligned_idxs & aligned_to_sexcat_idxs_rev
        non_mutual_hits = sexcat_to_aligned_idxs ^ aligned_to_sexcat_idxs_rev
        assert mutual_hits | non_mutual_hits == sexcat_to_aligned_idxs | aligned_to_sexcat_idxs_rev

        sexcat_in_non_mutual = set(i for i, j in non_mutual_hits)
        aligned_in_non_mutual = set(j for i, j in non_mutual_hits)
        exclusive_hits = set((i, j) for i, j in mutual_hits if i not in
                             sexcat_in_non_mutual and j not in aligned_in_non_mutual)

        print 'Non-mutual hits:', len(non_mutual_hits)
        print 'Mutual hits:', len(mutual_hits)
        print 'Exclusive hits:', len(exclusive_hits)

        self.non_mutual_hits = non_mutual_hits
        self.mutual_hits = mutual_hits
        self.exclusive_hits = exclusive_hits

    def extract_intensity_and_sequence_from_fastq(self,
                                                  fastq_fpath,
                                                  tile_key,
                                                  tile_num,
                                                  out_fpath,
                                                  hit_type='exclusive'):
        hit_list = list(getattr(self, hit_type + '_hits'))
        hit_aligned_idxs = [j for _, j in hit_list]
        rcs_coord_tups = set((tile_num, pt[0], pt[1]) for pt in self.rcs_in_frame)
        hit_given_rcs_coord_tup = \
                {(tile_num, pt[0], pt[1]): hit_list[i] for i, pt in enumerate(self.rcs_in_frame)}

        def flux_info_given_rcs_coord_tup(coord_tup):
            i, _ = hit_given_rcs_coord_tup[coord_tup]
            sexcat_pt = self.sexcat.points[i]
            return sexcat_pt.flux, sexcat_pt.flux_err

        with open(out_fpath, 'w') as out:
            for record in SeqIO.parse(open(fastq_fpath), 'fastq'):
                coord_tup = tuple(map(int, str(record.id).split(':')[-3:]))  # tile:r:c
                if coord_tup in rcs_coord_tups:
                    flux, flux_err = flux_info_given_rcs_coord_tup(coord_tup)
                    #out.write('\t'.join([record.id, self.image_data.fname,
                    record.description += ' Flux:%f Flux_err:%f' % (flux, flux_err)
                    SeqIO.write(record, out, 'fastq')
