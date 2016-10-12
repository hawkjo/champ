import os
import glob
import re
import h5py
import misc
from champ import hdf5tools
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import logging

log = logging.getLogger(__name__)


class IntensityScores(object):
    def __init__(self, h5_fpaths):
        """Initialize h5_fpaths and scores. scores is a dict accessed as:

            scores[h5_fpath][channel][pos_tup][read_name]
        """
        self.h5_fpaths = h5_fpaths
        self.raw_scores = {
            h5_fpath: {channel: {} for channel in hdf5tools.load_channel_names(h5_fpath)}
            for h5_fpath in h5_fpaths
            }
        self.scores = self.raw_scores

    def _make_isimportant_function(self, important_read_names):
        # Set cluster skip test
        if important_read_names == 'all':
            def isimportant(*args):
                return True
        else:
            if not isinstance(important_read_names, set):
                important_read_names = set(important_read_names)

            def isimportant(read_name):
                return read_name in important_read_names
        return isimportant

    def get_LDA_scores(self,
                       results_dirs,
                       lda_weights_fpath,
                       side_px=3,
                       verbose=True,
                       important_read_names='all'):
        isimportant = self._make_isimportant_function(important_read_names)

        # Read scores
        lda_weights = np.loadtxt(lda_weights_fpath)
        im_loc_re = re.compile('Channel_(.+)_Pos_(\d+)_(\d+)_')
        image_parsing_regex = re.compile(r'^(?P<channel>.+)_(?P<minor>\d+)_(?P<major>\d+)_')
        for h5_fpath, results_dir in zip(self.h5_fpaths, results_dirs):
            results_fpaths = glob.glob(os.path.join(results_dir, '*_all_read_rcs.txt'))
            if verbose:
                print h5_fpath
                print 'Num results files:', len(results_fpaths)

            for i, rfpath in enumerate(results_fpaths):
                rfname = os.path.basename(rfpath)
                try:
                    m = im_loc_re.match(rfname)
                    channel = m.group(1)
                    minor, major = tuple(int(m.group(i)) for i in (2, 3))
                except:
                    try:
                        m = image_parsing_regex.match(rfname)
                        channel = m.group('channel')
                        minor, major = int(m.group('minor')), int(m.group('major'))
                    except:
                        print rfname
                        raise

                pos_key = hdf5tools.get_image_key(major, minor)

                self.scores[h5_fpath][channel][(major, minor)] = {}

                with h5py.File(h5_fpath) as f:
                    im = np.array(f[channel][pos_key])

                for line in open(rfpath):
                    read_name, r, c = line.strip().split()
                    if not isimportant(read_name):
                        continue
                    r, c = map(misc.stoftoi, (r, c))
                    if (side_px <= r < im.shape[0] - side_px - 1
                        and side_px <= c < im.shape[0] - side_px - 1):
                        x = im[r - side_px:r + side_px + 1, c - side_px:c + side_px + 1].astype(np.float)
                        score = float(np.multiply(lda_weights, x).sum())
                        self.scores[h5_fpath][channel][(major, minor)][read_name] = score

    def normalize_scores(self, verbose=True):
        """Normalizes scores. The normalizing constant for each image is determined by

            Z = mode(pixel values) / median(all modes in h5_fpath)
        """

        def get_mode_in_im(im):
            w = 200
            hw = w / 2
            rmid, cmid = int(im.shape[0] / 2), int(im.shape[1] / 2)
            vmin, vmax = im.min(), im.max()
            # remove saturation
            pct95 = vmin + 0.95 * (vmax - vmin)
            vals = [v for v in im[rmid - hw:rmid + hw, cmid - hw:cmid + hw].flatten() if v < pct95]
            return misc.get_mode(vals)

        self.scores = {
            h5_fpath: {channel: {} for channel in hdf5tools.load_channel_names(h5_fpath)}
            for h5_fpath in self.h5_fpaths
            }
        self.normalizing_constants = {
            h5_fpath: {channel: {} for channel in hdf5tools.load_channel_names(h5_fpath)}
            for h5_fpath in self.h5_fpaths
            }
        for h5_fpath in self.h5_fpaths:
            if verbose: print os.path.basename(h5_fpath)
            for channel in self.scores[h5_fpath].keys():
                mode_given_pos_tup = {}
                for pos_tup in self.raw_scores[h5_fpath][channel].keys():
                    pos_key = hdf5tools.get_image_key(*pos_tup)
                    with h5py.File(h5_fpath) as f:
                        im = np.array(f[channel][pos_key])

                    mode_given_pos_tup[pos_tup] = get_mode_in_im(im)

                median_of_modes = np.median(mode_given_pos_tup.values())
                for pos_tup in mode_given_pos_tup.keys():
                    Z = mode_given_pos_tup[pos_tup] / float(median_of_modes)
                    self.normalizing_constants[h5_fpath][channel][pos_tup] = Z
                    im_scores = self.raw_scores[h5_fpath][channel][pos_tup]
                    self.scores[h5_fpath][channel][pos_tup] = {
                        read_name: im_scores[read_name] / Z
                        for read_name in self.get_read_names_in_image(h5_fpath, channel, pos_tup)
                        }
            if verbose: print

    def normalize_scores_by_ref_read_names(self, ref_read_names_given_channel, verbose=True):
        """Normalizes scores. The normalizing constant for each image is determined by

            Z = median(reference read scores) / 100
        """
        self.scores = {
            h5_fpath: {channel: {} for channel in hdf5tools.load_channel_names(h5_fpath)}
            for h5_fpath in self.h5_fpaths
            }
        self.normalizing_constants = {
            h5_fpath: {channel: {} for channel in hdf5tools.load_channel_names(h5_fpath)}
            for h5_fpath in self.h5_fpaths
            }
        for h5_fpath in self.h5_fpaths:
            log.debug(os.path.basename(h5_fpath))
            for channel in self.scores[h5_fpath].keys():
                ref_read_names = ref_read_names_given_channel[channel]
                for pos_tup in self.raw_scores[h5_fpath][channel].keys():
                    ref_read_names_in_image = (
                        self.get_read_names_in_image(h5_fpath, channel, pos_tup)
                        & ref_read_names
                    )
                    if len(ref_read_names_in_image) < 10:
                        print 'Warning: 10 > {} reference reads in im_idx {}'.format(
                            len(ref_read_names_in_image), (h5_fpath, channel, pos_tup)
                        )

                    med = np.median(
                        [self.raw_scores[h5_fpath][channel][pos_tup][read_name]
                         for read_name in ref_read_names_in_image]
                    )

                    Z = med / 100.0
                    self.normalizing_constants[h5_fpath][channel][pos_tup] = Z
                    im_scores = self.raw_scores[h5_fpath][channel][pos_tup]
                    self.scores[h5_fpath][channel][pos_tup] = {
                        read_name: im_scores[read_name] / Z
                        for read_name in self.get_read_names_in_image(h5_fpath, channel, pos_tup)
                        }

    def get_read_names_in_image(self, h5_fpath, channel, pos_tup):
        return set(self.raw_scores[h5_fpath][channel][pos_tup].keys())

    def build_score_given_read_name_given_channel(self):
        self.score_given_read_name_in_channel = {
            h5_fpath: {channel: {} for channel in hdf5tools.load_channel_names(h5_fpath)}
            for h5_fpath in self.h5_fpaths
            }
        for h5_fpath in self.h5_fpaths:
            print h5_fpath
            i = 0
            for channel in self.scores[h5_fpath].keys():
                score_given_read_name = self.score_given_read_name_in_channel[h5_fpath][channel]
                for pos_tup in self.scores[h5_fpath][channel].keys():
                    for read_name, score in self.scores[h5_fpath][channel][pos_tup].items():
                        score_given_read_name[read_name] = score
                        i += 1

    def plot_normalization_constants(self):
        for h5_fpath in self.h5_fpaths:
            nMajor_pos, nminor_pos = hdf5tools.calculate_grid_dimensions(h5_fpath)
            for channel in sorted(self.scores[h5_fpath].keys()):
                fig, ax = plt.subplots(figsize=(10, 1))
                M = np.empty((nminor_pos + 1, nMajor_pos + 1))
                M[:] = None

                for pos_tup in self.scores[h5_fpath][channel].keys():
                    col, row = pos_tup
                    M[row, col] = self.normalizing_constants[h5_fpath][channel][pos_tup]

                ms = ax.matshow(M)
                cbar = plt.colorbar(ms)

                ax.set_title('Normalizing constants in {} Channel {}'.format(os.path.basename(h5_fpath), channel))
                ax.set_aspect(1)
                ax.xaxis.set_ticks_position('bottom')

    def plot_aligned_images(self, colors='rgbcmyk', markers='o*^sv+x'):
        n = len(self.h5_fpaths)
        fig, axes = plt.subplots(n, 1, figsize=(10, n))

        for h5_fpath, ax in zip(self.h5_fpaths, axes):
            for channel, color, marker in zip(
                    sorted(self.scores[h5_fpath].keys()), colors, markers
            ):
                rs, cs = [], []
                for pos_tup in self.scores[h5_fpath][channel].keys():
                    c, r = pos_tup
                    rs.append(r)
                    cs.append(c)
                ax.plot(cs, rs, marker, color=color, alpha=0.4, label=channel)

            ax.set_title('Aligned images in {}'.format(os.path.basename(h5_fpath)))
            nMajor_pos, nminor_pos = hdf5tools.calculate_grid_dimensions(h5_fpath)
            ax.set_ylim((nminor_pos, -1))
            ax.set_xlim((-1, 1.15 * nMajor_pos))  # Add room for legend
            ax.set_aspect(1)
            ax.legend()
        return fig, axes

    def print_reads_per_channel(self):
        reads_in_channel = defaultdict(set)
        for h5_fpath in self.h5_fpaths:
            for channel in self.scores[h5_fpath].keys():
                for score_given_read_name in self.scores[h5_fpath][channel].values():
                    rnames = score_given_read_name.keys()
                    print("len(read_names) = %d" % len(rnames))
                    reads_in_channel[channel].update(rnames)
        for channel, read_names in sorted(reads_in_channel.items()):
            print 'All reads found in channel {}: {:,d}'.format(channel, len(read_names))

    def build_good_read_names(self, good_num_ims_cutoff, pos_tups_cutoff=1):
        pos_tups_given_read_name = defaultdict(set)
        h5_fpaths_given_read_name = defaultdict(set)
        for h5_fpath in self.h5_fpaths:
            for channel in self.scores[h5_fpath].keys():
                for pos_tup in self.scores[h5_fpath][channel].keys():
                    for read_name in self.scores[h5_fpath][channel][pos_tup].keys():
                        pos_tups_given_read_name[read_name].add(pos_tup)
                        h5_fpaths_given_read_name[read_name].add(h5_fpath)
        self.good_read_names = set(
            read_name for read_name, pos_tups in pos_tups_given_read_name.items()
            if len(h5_fpaths_given_read_name[read_name]) >= good_num_ims_cutoff
        )

    def write_values_by_seq(self,
                            course_trait_name,
                            course_trait_list,
                            h5_fpaths,
                            attrs_dict,
                            seqs_of_interest,
                            read_names_given_seq,
                            channel_of_interest,
                            out_fpath,
                            ):
        """
        Writes output in array-like format.

        Params:
            :str:   course_trait_name - description of defining trait for h5_fpaths
            :list:  course_trait_list - list of defining trait values for each h5_fpath
            :list:  h5_fpaths - the desired subset of IntensityScores.h5_fpaths to output
            :dict:  attrs_dict - a dict containing additional (str) attributes of interest
            :key:   channel of interest - key for channel of interest
            :iter:  seqs_of_interest
            :dict:  read_names_given_seq
            :str:   out_fpath
        """
        assert len(h5_fpaths) == len(course_trait_list), (h5_fpaths, course_trait_list)
        with open(out_fpath, 'w') as out:
            out.write('# Defining Course Trait: {}\n'.format(course_trait_name))
            out.write('\t'.join(map(str, course_trait_list)) + '\n')
            out.write('# HDF5 Files\n')
            out.write('\n'.join(h5_fpaths) + '\n')
            out.write('# Channel: {}\n'.format(str(channel_of_interest)))
            for k, v in sorted(attrs_dict.items()):
                out.write('# {}: {}\n'.format(k, v))
            for seq in seqs_of_interest:
                out.write(seq + '\n')
                read_names = list(read_names_given_seq[seq])
                out.write('\t'.join(read_names) + '\n')
                for h5_fpath in h5_fpaths:
                    score_given_read_name = self.score_given_read_name_in_channel[h5_fpath][channel_of_interest]
                    out.write('\t'.join(str(float(score_given_read_name[read_name]))
                                        if read_name in score_given_read_name
                                        else '-'
                                        for read_name in read_names)
                              + '\n')
