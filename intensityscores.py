import os
import glob
import string
import nd2tools
import misc
import misctools
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict, Counter

class IntensityScores(object):
    def __init__(self, nd2s):
        """Initialize nd2s and scores. scores is a dict accessed as:

            scores[nd2][im_idx][read_name]
        """
        self.nd2s = nd2s
        self.raw_scores = {nd2: {} for nd2 in nd2s}
        self.scores = self.raw_scores

    def get_LDA_scores(self,
                       results_dir_given_nd2,
                       lda_weights_fpath,
                       side_px=3,
                       verbose=True,
                       important_read_names='all'):

        # Set cluster skip test
        if important_read_names == 'all':
            def isimportant(*args):
                return True
        else:
            if not isinstance(important_read_names, set):
                important_read_names = set(important_read_names)
            def isimportant(read_name):
                return read_name in important_read_names

        # Read scores
        lda_weights = np.loadtxt(lda_weights_fpath)
        for nd2 in self.nd2s:
            results_dir = results_dir_given_nd2[nd2]
            assert os.path.isdir(results_dir), results_dir
            results_fpaths = glob.glob(os.path.join(results_dir, '*_all_read_rcs.txt'))
            if verbose: print 'Num results files:', len(results_fpaths)
        
            for i, rfpath in enumerate(results_fpaths):
                rfname = os.path.basename(rfpath)
                im_idx = int(rfname[:rfname.index('_')])
                self.scores[nd2][im_idx] = {}
                if verbose: misctools.dot()
        
                im = nd2[im_idx].data
                for line in open(rfpath):
                    read_name, r, c = line.strip().split()
                    if not isimportant(read_name):
                        continue
                    r, c = map(misc.stoftoi, (r, c))
                    if (side_px <= r < im.shape[0] - side_px - 1
                        and side_px <= c < im.shape[0] - side_px - 1):
                        x = im[r-side_px:r+side_px+1, c-side_px:c+side_px+1].astype(np.float)
                        score = np.multiply(lda_weights, x).sum()
                        self.scores[nd2][im_idx][read_name] = score
            if verbose: print

    def normalize_scores(self, reference_read_names, verbose=True):
        """Normalizes scores. The normalizing constant for each image is determined by
            
            Z = median(reference read scores bounded below) / median(all medians in nd2)

        where 'bounded below' means read scores are artificially set to 1 if otherwise lower.
        """
        if not isinstance(reference_read_names, set):
            reference_read_names = set(reference_read_names)

        self.scores = {nd2: {} for nd2 in self.nd2s}
        self.normalizing_constants = {nd2: {} for nd2 in self.nd2s}
        for nd2 in self.nd2s:
            if verbose: print nd2tools.bname_given_nd2(nd2)
            median_given_im_idx = {}
            for im_idx in self.raw_scores[nd2].keys():
                if verbose: misctools.dot()
                reference_read_names_in_image = (self.get_read_names_in_image(nd2, im_idx)
                                                 & reference_read_names)
                median_given_im_idx[im_idx]  = np.median(
                    [max(self.raw_scores[nd2][im_idx][read_name], 1)
                     for read_name in reference_read_names_in_image]
                )

            median_of_medians = np.median(median_given_im_idx.values())
            for im_idx in self.raw_scores[nd2].keys():
                Z = median_given_im_idx[im_idx] / float(median_of_medians)
                self.normalizing_constants[nd2][im_idx] = Z
                im_scores = self.raw_scores[nd2][im_idx]
                self.scores[nd2][im_idx] = {
                    read_name: im_scores[read_name] / Z
                    for read_name in self.get_read_names_in_image(nd2, im_idx)
                }
        
    def get_read_names_in_image(self, nd2, im_idx):
        return set(self.raw_scores[nd2][im_idx].keys())

    def plot_normalization_constants(self):
        for nd2 in self.nd2s:
            def spawn_matrix():
                M = np.empty(nd2tools.nrows_and_ncols(nd2))
                M[:] = None
                return M
            M_given_channel_idx = defaultdict(spawn_matrix)

            for im_idx in self.scores[nd2].keys():
                channel_idx = im_idx % len(nd2.channels)
                pos_name = nd2tools.convert_nd2_coordinates(nd2, im_idx=im_idx, outfmt='pos_name')
                row = string.uppercase.index(pos_name[0])
                col = int(pos_name[1:])
                M_given_channel_idx[channel_idx][row, col] = self.normalizing_constants[nd2][im_idx]

            for channel_idx, M in M_given_channel_idx.items():
                fig, ax = plt.subplots(figsize=(20, 1))
                ms = ax.matshow(M)
                cbar = plt.colorbar(ms)

                ax.set_title('Normalizing constants in %s channel %s'
                             % (nd2tools.bname_given_nd2(nd2), nd2.channels[channel_idx]))
                ax.set_aspect(1)
                ax.set_xlim(ax.get_xlim()[::-1])
                ax.set_ylim(ax.get_ylim()[::-1])
                ax.set_yticks(range(M.shape[0]))
                ax.set_yticklabels(string.uppercase[:M.shape[0]])
                ax.xaxis.set_ticks_position('bottom')

    def plot_aligned_images(self):
        for nd2 in self.nd2s:
            rs = defaultdict(list)
            cs = defaultdict(list)
            for im_idx in self.scores[nd2].keys():
                channel_idx = im_idx % len(nd2.channels)
                r, c = nd2tools.convert_nd2_coordinates(nd2, im_idx=im_idx, outfmt='pos_coords')
                rs[channel_idx].append(r)
                cs[channel_idx].append(c)

            fig, ax = plt.subplots(figsize=(10, 7))
            for channel_idx, marker in zip(rs.keys(), 'o*x^'):
                color = nd2.channels[channel_idx].lower()
                ax.plot(cs[channel_idx], rs[channel_idx], marker,
                        color=color, alpha=0.4,
                        label=nd2.channels[channel_idx])
            ax.set_title('Aligned images in %s' % nd2tools.bname_given_nd2(nd2))
            ax.set_aspect(1)
            ax.set_ylim(ax.get_ylim()[::-1])
            xlim = ax.get_xlim()
            ax.set_xlim((xlim[0], xlim[1] + 0.15 * (xlim[1] - xlim[0])))  # Add room for legend
            ax.legend()

    def print_reads_per_channel(self):
        reads_per_channel = Counter()
        for nd2 in self.nd2s:
            for im_idx, score_given_read_name in self.scores[nd2].items():
                channel_idx = im_idx % len(nd2.channels)
                reads_per_channel[channel_idx] += len(score_given_read_name)
        for channel_idx, num_reads in sorted(reads_per_channel.items()):
            print 'All reads found in channel %d: %d' % (channel_idx, num_reads)

    def build_good_read_names(self, good_num_ims_cutoff):
        pos_names_given_read_name = defaultdict(set)
        num_nd2s_given_read_name = defaultdict(set)
        for nd2 in self.nd2s:
            for im_idx in self.scores[nd2].keys():
                pos_name = nd2tools.convert_nd2_coordinates(nd2, outfmt='pos_name', im_idx=im_idx)
                for read_name in self.scores[nd2][im_idx].keys():
                    pos_names_given_read_name[read_name].add(pos_name)
                    num_nd2s_given_read_name[read_name].add(nd2)
        self.good_read_names = set(
            read_name for read_name, pos_names in pos_names_given_read_name.items()
            if len(pos_names) == 1
            and len(num_nd2s_given_read_name[read_name]) >= good_num_ims_cutoff
        )
