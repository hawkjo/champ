import os
import glob
import gc
import re
import h5py
import misc
from champ import hdf5tools
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import biofits
import multiprocessing
import logging
from multiprocessing import Process
from multiprocessing.queues import SimpleQueue

log = logging.getLogger(__name__)


def get_image_mode(image):
    w = 200
    hw = w / 2
    rmid, cmid = int(image.shape[0] / 2), int(image.shape[1] / 2)
    vmin, vmax = image.min(), image.max()
    # remove saturation
    pct95 = vmin + 0.95 * (vmax - vmin)
    vals = [v for v in image[rmid - hw:rmid + hw, cmid - hw:cmid + hw].flatten() if v < pct95]
    return misc.get_mode(vals)


def _thread_get_normalization_constants(h5_filename, results_queue):
    normalization_constant = {}
    with h5py.File(h5_filename, 'r') as h5:
        for channel, fov_and_image in h5.items():
            normalization_constant[channel] = {}
            image_modes = {}
            for fov, image in fov_and_image.items():
                image_modes[fov] = get_image_mode(image.value)
            median_of_modes = np.median(image_modes.values())
            for fov, mode in image_modes.items():
                normalization_constant[channel][fov] = mode / median_of_modes
    results_queue.put((h5_filename, normalization_constant))


def get_normalization_constants(h5_filenames):
    """
    Determines per-concentration pixel intensity normalization factors for each field of view.
    Does not perform normalization!

    """
    results_queue = SimpleQueue()
    processes = []
    for h5_filename in h5_filenames:
        p = Process(target=_thread_get_normalization_constants, args=(h5_filename, results_queue))
        processes.append(p)
        p.start()
    normalization_constant = {}
    for _ in h5_filenames:
        h5_filename, constants = results_queue.get()
        normalization_constant[h5_filename] = constants
    for p in processes:
        p.join()
    return normalization_constant


def _thread_calculate_lda_scores(h5_fpath, results_dir, normalization_constants, lda_weights, results_queue):
    side_pixels = 3
    image_parsing_regex = re.compile(r'^(?P<channel>.+)_(?P<minor>\d+)_(?P<major>\d+)_')
    scores = defaultdict(dict)
    results_paths = glob.glob(os.path.join(results_dir, '*_all_read_rcs.txt'))
    for i, rfpath in enumerate(results_paths):
        rfname = os.path.basename(rfpath)
        m = image_parsing_regex.match(rfname)
        channel = m.group('channel')
        if channel not in scores:
            scores[channel] = defaultdict(list)
        minor, major = int(m.group('minor')), int(m.group('major'))
        field_of_view = hdf5tools.get_image_key(major, minor)
        norm_constant = normalization_constants[channel][field_of_view]
        with h5py.File(h5_fpath) as f:
            im = f[channel][field_of_view].value

        with open(rfpath) as f:
            for line in f:
                read_name, r, c = line.strip().split()
                r, c = map(misc.stoftoi, (r, c))
                far_enough_from_vertical_edges = side_pixels <= r < im.shape[0] - side_pixels - 1
                far_enough_from_horizontal_edges = side_pixels <= c < im.shape[0] - side_pixels - 1
                if not far_enough_from_vertical_edges or not far_enough_from_horizontal_edges:
                    continue
                x = im[r - side_pixels:r + side_pixels + 1, c - side_pixels:c + side_pixels + 1].astype(np.float)
                score = float(np.multiply(lda_weights, x).sum()) * norm_constant
                # If a read shows up in multiple fields of view, we'll want to just take the mean value
                scores[channel][read_name].append(score)

    # all our data is in lists, with most lists having only a single member
    # if a particular read was observed in more than one field of view, we take the mean of its values
    mean_scores = {}
    for channel, read_name_data in scores.items():
        mean_scores[channel] = {}
        for read_name, scores in read_name_data.items():
            mean_score = np.mean(scores) if len(scores) > 1 else scores[0]
            mean_scores[channel][read_name] = mean_score
    results_queue.put((h5_fpath, mean_scores))


class LDAScores(object):
    """
    This object simply marshals intensity data for a single DNA cluster from a multi-tier dictionary. Each DNA cluster
    may or may not have an intensity value for each concentration under each color. We handle the fact that gaps exist
    here by slotting in data as it arrives. Later, when we fit intensity curves, we can figure out which concentrations
    were missing as they will have None instead of a float.

    """
    __slots__ = '_conditions', '_scores'

    def __init__(self, conditions, channels):
        self._conditions = conditions
        self._scores = {channel: [None for _ in conditions] for channel in channels}

    def add(self, channel, condition, score):
        index = self._conditions.index(condition)
        self._scores[channel][index] = score

    def get(self, channel):
        # it's not guaranteed that each read will appear in all channels
        return self._scores.get(channel)


def calculate_lda_scores(h5_paths, results_directories, normalization_constants, lda_weights_path):
    processes = []
    results_queue = SimpleQueue()
    lda_weights = np.loadtxt(lda_weights_path)
    for h5_fpath, results_dir in zip(h5_paths, results_directories):
        p = Process(target=_thread_calculate_lda_scores, args=(h5_fpath,
                                                               results_dir,
                                                               normalization_constants[h5_fpath],
                                                               lda_weights,
                                                               results_queue))
        processes.append(p)
        p.start()

    lda_scores = {}
    for _ in h5_paths:
        h5_filename, score_data = results_queue.get()
        channels = score_data.keys()
        for channel, read_name_scores in score_data.items():
            for read_name, score in read_name_scores.items():
                if read_name not in lda_scores:
                    lda_scores[read_name] = LDAScores(h5_paths, channels)
                lda_scores[read_name].add(channel, h5_filename, score)
        gc.collect()
    for p in processes:
        p.join()
    return lda_scores


def get_reasonable_process_count():
    # Determines how many processes to use for parallel tasks. For small machines, it uses all cores, but on
    # larger systems it will leave a few cores available for other processes under the assumption
    # that it needs to perform background functions and not interfere with SSH and such
    try:
        core_count = multiprocessing.cpu_count()
    except NotImplementedError:
        print("Could not determine number of cores. Using one process only!")
        return 1
    if core_count <= 4:
        return core_count
    if core_count <= 16:
        return core_count - 1
    return core_count - 4


def _thread_calculate_kds(concentrations, lda_scores, channel, results_queue):
    results = {}
    update_level = int(len(lda_scores) / 10)
    for n, (read_name, scores) in enumerate(lda_scores.items()):
        if n % update_level == 0:
            print("%.1f%% done." % (100.0 * float(n) / len(lda_scores)))
        read_concentrations = []
        read_intensities = []
        for concentration, score in zip(concentrations, scores.get(channel)):
            if score is not None:
                read_concentrations.append(concentration)
                read_intensities.append(score)
        if len(read_intensities) < 4:
            continue
        try:
            _, _, _, _, kd, kd_stddev = biofits.fit_hyperbola(read_concentrations, read_intensities)
            results[read_name] = (kd, kd_stddev)
        except RuntimeError:
            continue
    results_queue.put(results)


def calculate_kds(h5_paths, lda_scores, channel):
    import random
    concentrations = [misc.parse_concentration(fpath) for fpath in h5_paths]
    results_queue = SimpleQueue()
    process_count = get_reasonable_process_count()
    print("Using %d cores" % process_count)
    kds = {}
    print("Calculating %d KDs" % len(lda_scores))
    split_scores = [{} for _ in range(process_count)]
    for n, (read_name, lda_score) in enumerate(lda_scores.items()):
        if random.random() < 0.005:
            split_scores[n % process_count][read_name] = lda_score

    processes = []
    for split_score in split_scores:
        print("Sending %d LDA scores to a thread" % len(split_score))
        p = Process(target=_thread_calculate_kds, args=(concentrations, split_score, channel, results_queue))
        processes.append(p)
        p.start()
    print("Waiting for results")
    for _ in range(process_count):
        results = results_queue.get()
        kds.update(results)
    for p in processes:
        p.join()
    return kds


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
        if n == 1:
            # axes is a tuple unless n==1, but the rest of the code depends on it being iterable
            axes = (axes,)
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

    def print_reads_per_channel(self):
        reads_in_channel = defaultdict(set)
        for h5_fpath in self.h5_fpaths:
            for channel in self.scores[h5_fpath].keys():
                for score_given_read_name in self.scores[h5_fpath][channel].values():
                    reads_in_channel[channel].update(score_given_read_name.keys())
        for channel, read_names in sorted(reads_in_channel.items()):
            print 'All reads found in channel {}: {:,d}'.format(channel, len(read_names))

    def build_good_read_names(self, good_num_ims_cutoff):
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
