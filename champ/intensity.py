import glob
import logging
import os
import re
from collections import defaultdict
import h5py
import matplotlib.pyplot as plt
import numpy as np
from champ import grid
from champ import hdf5tools, misc
from scipy.optimize import minimize
from sklearn.neighbors import KernelDensity

log = logging.getLogger(__name__)


def calculate_intensities(h5_filenames):
    results = defaultdict(list)
    for filename in h5_filenames:
        h5 = h5py.File(filename, "r")
        g = grid.GridImages(h5, "NGS_blue")
        for image in g.left_iter():
            results[filename].append(np.mean(image))
    return results


def plot_timecourse(times, intensities):
    plt.boxplot(intensities)
    plt.xticks(range(len(intensities)), times)
    plt.show()


class IntensityScores(object):
    def __init__(self, h5_filepaths):
        """Initialize h5_filepaths and scores. scores is a dict accessed as:

            scores[h5_fpath][channel][pos_tup][read_name]
        """
        self.h5_filepaths = h5_filepaths
        self.raw_scores = {
            h5_fpath: {channel: {} for channel in hdf5tools.load_channel_names(h5_fpath)}
            for h5_fpath in h5_filepaths
            }
        self.scores = self.raw_scores

    def get_LDA_scores(self, results_dirs, lda_weights_fpath, side_px=3, important_read_names='all'):
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
        image_parsing_regex = re.compile(r'^(?P<channel>.+)_(?P<minor>\d+)_(?P<major>\d+)_')
        for h5_fpath, results_dir in zip(self.h5_filepaths, results_dirs):
            results_fpaths = glob.glob(os.path.join(results_dir, '*_all_read_rcs.txt'))
            print('Num results files: %d' % len(results_fpaths))

            for i, result_path in enumerate(results_fpaths):
                result_filename = os.path.basename(result_path)
                m = image_parsing_regex.match(result_filename)
                channel = m.group('channel')
                minor, major = int(m.group('minor')), int(m.group('major'))
                position = hdf5tools.get_image_key(major, minor)
                self.scores[h5_fpath][channel][(major, minor)] = {}

                with h5py.File(h5_fpath) as f:
                    im = np.array(f[channel][position])

                with open(result_path) as f:
                    for line in f:
                        read_name, row, column = line.strip().split()
                        if not isimportant(read_name):
                            continue
                        row, column = map(misc.stoftoi, (row, column))
                        if side_px <= row < im.shape[0] - side_px - 1 and side_px <= column < im.shape[0] - side_px - 1:
                            x = im[row - side_px:row + side_px + 1, column - side_px:column + side_px + 1].astype(np.float)
                            score = np.multiply(lda_weights, x).sum()
                            self.scores[h5_fpath][channel][(major, minor)][read_name] = score

    def get_mode(self, im):
        w = 200
        hw = w / 2
        rmid, cmid = int(im.shape[0] / 2), int(im.shape[1] / 2)
        vmin, vmax = im.min(), im.max()
        bandwidth = (vmax - vmin) / 200
        kdf = KernelDensity(bandwidth=bandwidth)
        # remove saturation
        pct95 = vmin + 0.95 * (vmax - vmin)
        vals = [v for v in im[rmid - hw:rmid + hw, cmid - hw:cmid + hw].flatten() if v < pct95]
        kdf.fit(np.array(vals).reshape(len(vals), 1))

        def neg_kdf(x):
            return -kdf.score(x)

        res = minimize(neg_kdf, x0=np.median(im.flatten()), method='Nelder-Mead')
        assert res.success, res
        return res.x

    def normalize_scores(self):
        """Normalizes scores. The normalizing constant for each image is determined by

            Z = median(reference read scores bounded below) / median(all medians in h5_fpath)

        where 'bounded below' means read scores are artificially set to 1 if otherwise lower.
        """

        self.scores = {
            h5_fpath: {channel: {} for channel in hdf5tools.load_channel_names(h5_fpath)}
            for h5_fpath in self.h5_filepaths
            }
        self.normalizing_constants = {
            h5_fpath: {channel: {} for channel in hdf5tools.load_channel_names(h5_fpath)}
            for h5_fpath in self.h5_filepaths
            }
        for h5_fpath in self.h5_filepaths:
            for channel in self.scores[h5_fpath].keys():
                mode_given_pos_tup = {}
                for pos_tup in self.raw_scores[h5_fpath][channel].keys():
                    pos_key = hdf5tools.get_image_key(*pos_tup)
                    with h5py.File(h5_fpath) as f:
                        im = np.array(f[channel][pos_key])
                    mode_given_pos_tup[pos_tup] = self.get_mode(im)

                median_of_modes = np.median(mode_given_pos_tup.values())
                for pos_tup in mode_given_pos_tup.keys():
                    Z = mode_given_pos_tup[pos_tup] / float(median_of_modes)
                    self.normalizing_constants[h5_fpath][channel][pos_tup] = Z
                    im_scores = self.raw_scores[h5_fpath][channel][pos_tup]
                    self.scores[h5_fpath][channel][pos_tup] = {
                        read_name: im_scores[read_name] / Z
                        for read_name in self.get_read_names_in_image(h5_fpath, channel, pos_tup)
                        }

    def normalize_scores_by_ref_read_names(self, ref_read_names_given_channel):
        """Normalizes scores. The normalizing constant for each image is determined by

            Z = median(reference read scores) / 100
        """
        self.scores = {
            h5_fpath: {channel: {} for channel in hdf5tools.load_channel_names(h5_fpath)}
            for h5_fpath in self.h5_filepaths
            }
        self.normalizing_constants = {
            h5_fpath: {channel: {} for channel in hdf5tools.load_channel_names(h5_fpath)}
            for h5_fpath in self.h5_filepaths
            }
        for h5_fpath in self.h5_filepaths:
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
            for h5_fpath in self.h5_filepaths
            }
        for h5_fpath in self.h5_filepaths:
            for channel in self.scores[h5_fpath].keys():
                score_given_read_name = self.score_given_read_name_in_channel[h5_fpath][channel]
                for pos_tup in self.scores[h5_fpath][channel].keys():
                    for read_name, score in self.scores[h5_fpath][channel][pos_tup].items():
                        score_given_read_name[read_name] = score

    def plot_normalization_constants(self):
        for h5_fpath in self.h5_filepaths:
            basename = os.path.basename(h5_fpath)
            num_columns, num_rows = hdf5tools.calculate_grid_dimensions(h5_fpath)
            for channel in sorted(self.scores[h5_fpath].keys()):
                M = np.empty((num_rows, num_columns))
                M[:] = None

                for pos_tup in self.scores[h5_fpath][channel].keys():
                    col, row = pos_tup
                    M[row, col] = self.normalizing_constants[h5_fpath][channel][pos_tup]

                fig, ax = plt.subplots(figsize=(20, 1))
                ms = ax.matshow(M)
                plt.colorbar(ms)

                ax.set_title('Normalizing constants in {} Channel {}'.format(basename, channel))
                ax.set_aspect(1)
                ax.xaxis.set_ticks_position('bottom')
                yield basename, channel, fig

    def plot_aligned_images(self, colors, markers):
        for h5_fpath in self.h5_filepaths:
            basename = os.path.basename(h5_fpath)
            fig, ax = plt.subplots(figsize=(10, 7))
            for channel, color, marker in zip(
                    sorted(self.scores[h5_fpath].keys()), colors, markers
            ):
                rs, cs = [], []
                for c, r in self.scores[h5_fpath][channel].keys():
                    rs.append(r)
                    cs.append(c)
                ax.plot(cs, rs, marker, color=color, alpha=0.4, label=channel)

            ax.set_title('Aligned images in {}'.format(basename))
            nMajor_pos, nminor_pos = hdf5tools.calculate_grid_dimensions(h5_fpath)
            ax.set_ylim((nminor_pos, -1))
            ax.set_xlim((-1, 1.15 * nMajor_pos))  # Add room for legend
            ax.set_aspect(1)
            ax.legend()
            yield basename, fig

    def print_reads_per_channel(self):
        reads_in_channel = defaultdict(set)
        for h5_fpath in self.h5_filepaths:
            for channel in self.scores[h5_fpath].keys():
                for score_given_read_name in self.scores[h5_fpath][channel].values():
                    reads_in_channel[channel].update(score_given_read_name.keys())
        for channel, read_names in sorted(reads_in_channel.items()):
            print('All reads found in channel {}: {:,d}'.format(channel, len(read_names)))

    def build_good_read_names(self, good_num_ims_cutoff):
        pos_tups_given_read_name = defaultdict(set)
        h5_filepaths_given_read_name = defaultdict(set)
        for h5_fpath in self.h5_filepaths:
            for channel in self.scores[h5_fpath].keys():
                for pos_tup in self.scores[h5_fpath][channel].keys():
                    for read_name in self.scores[h5_fpath][channel][pos_tup].keys():
                        pos_tups_given_read_name[read_name].add(pos_tup)
                        h5_filepaths_given_read_name[read_name].add(h5_fpath)
        self.good_read_names = set(
            read_name for read_name, pos_names in pos_tups_given_read_name.items()
            if len(pos_names) == 1
            and len(h5_filepaths_given_read_name[read_name]) >= good_num_ims_cutoff
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
