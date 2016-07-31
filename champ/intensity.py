import re
import h5py
from champ import hdf5tools, constants, error, projectinfo
from sklearn.neighbors import KernelDensity
import logging
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import random
import itertools
import misc
from collections import defaultdict
from scipy.optimize import minimize, curve_fit
from matplotlib.ticker import MultipleLocator
import functools


log = logging.getLogger(__name__)


class IntensityScores(object):
    def __init__(self, h5_filepaths):
        """Initialize h5_filepaths and scores. scores is a dict accessed as:

            scores[h5_fpath][channel][pos_tup][read_name]
        """
        self.h5_filepaths = h5_filepaths
        self.raw_scores = {
            h5_fpath: {channel: {} for channel in hdf5tools.get_channel_names(h5_fpath)}
            for h5_fpath in h5_filepaths
            }
        self.scores = self.raw_scores

    def get_LDA_scores(self,
                       results_dirs,
                       lda_weights_fpath,
                       side_px=3,
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
        im_loc_re = re.compile('Channel_(.+)_Pos_(\d+)_(\d+)_')
        for h5_fpath, results_dir in zip(self.h5_filepaths, results_dirs):
            results_fpaths = glob.glob(os.path.join(results_dir, '*_all_read_rcs.txt'))
            log.debug('Num results files: %d' % len(results_fpaths))

            for i, rfpath in enumerate(results_fpaths):
                rfname = os.path.basename(rfpath)
                m = im_loc_re.match(rfname)
                channel = m.group(1)
                pos_tup = tuple(int(m.group(i)) for i in (2, 3))
                pos_key = hdf5tools.dset_name_given_coords(*pos_tup)
                self.scores[h5_fpath][channel][pos_tup] = {}

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
                        score = np.multiply(lda_weights, x).sum()
                        self.scores[h5_fpath][channel][pos_tup][read_name] = score

    def normalize_scores(self):
        """Normalizes scores. The normalizing constant for each image is determined by

            Z = median(reference read scores bounded below) / median(all medians in h5_fpath)

        where 'bounded below' means read scores are artificially set to 1 if otherwise lower.
        """

        def get_mode(im):
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

        self.scores = {
            h5_fpath: {channel: {} for channel in hdf5tools.get_channel_names(h5_fpath)}
            for h5_fpath in self.h5_filepaths
            }
        self.normalizing_constants = {
            h5_fpath: {channel: {} for channel in hdf5tools.get_channel_names(h5_fpath)}
            for h5_fpath in self.h5_filepaths
            }
        for h5_fpath in self.h5_filepaths:
            log.debug(os.path.basename(h5_fpath))
            for channel in self.scores[h5_fpath].keys():
                mode_given_pos_tup = {}
                for pos_tup in self.raw_scores[h5_fpath][channel].keys():
                    pos_key = hdf5tools.dset_name_given_coords(*pos_tup)
                    with h5py.File(h5_fpath) as f:
                        im = np.array(f[channel][pos_key])

                    mode_given_pos_tup[pos_tup] = get_mode(im)

                median_of_modes = np.median(mode_given_pos_tup.values())
                for pos_tup in mode_given_pos_tup.keys():
                    Z = mode_given_pos_tup[pos_tup] / float(median_of_modes)
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
            h5_fpath: {channel: {} for channel in hdf5tools.get_channel_names(h5_fpath)}
            for h5_fpath in self.h5_filepaths
            }
        for h5_fpath in self.h5_filepaths:
            log.debug(h5_fpath)
            for channel in self.scores[h5_fpath].keys():
                score_given_read_name = self.score_given_read_name_in_channel[h5_fpath][channel]
                for pos_tup in self.scores[h5_fpath][channel].keys():
                    for read_name, score in self.scores[h5_fpath][channel][pos_tup].items():
                        score_given_read_name[read_name] = score

    def plot_normalization_constants(self):
        for h5_fpath in self.h5_filepaths:
            nMajor_pos, nminor_pos = hdf5tools.get_nMajor_nminor_pos(h5_fpath)
            for channel in sorted(self.scores[h5_fpath].keys()):
                M = np.empty((nminor_pos, nMajor_pos))
                M[:] = None

                for pos_tup in self.scores[h5_fpath][channel].keys():
                    col, row = pos_tup
                    M[row, col] = self.normalizing_constants[h5_fpath][channel][pos_tup]

                fig, ax = plt.subplots(figsize=(20, 1))
                ms = ax.matshow(M)
                plt.colorbar(ms)

                ax.set_title('Normalizing constants in {} Channel {}'.format(os.path.basename(h5_fpath), channel))
                ax.set_aspect(1)
                ax.xaxis.set_ticks_position('bottom')

    def plot_aligned_images(self, colors, markers):
        for h5_fpath in self.h5_filepaths:
            fig, ax = plt.subplots(figsize=(10, 7))
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
            nMajor_pos, nminor_pos = hdf5tools.get_nMajor_nminor_pos(h5_fpath)
            ax.set_ylim((nminor_pos, -1))
            ax.set_xlim((-1, 1.15 * nMajor_pos))  # Add room for legend
            ax.set_aspect(1)
            ax.legend()

    def print_reads_per_channel(self):
        reads_in_channel = defaultdict(set)
        for h5_fpath in self.h5_filepaths:
            for channel in self.scores[h5_fpath].keys():
                for score_given_read_name in self.scores[h5_fpath][channel].values():
                    reads_in_channel[channel].update(score_given_read_name.keys())
        for channel, read_names in sorted(reads_in_channel.items()):
            log.debug('All reads found in channel {}: {:,d}'.format(channel, len(read_names)))

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


def plot_single_mismatch_ddgs(seq_ddGs, seq_ddG_error, target, reference_sequence_name, fs=16):
    base_colors = dict(A='b', C='darkgoldenrod', G='g', T='r')
    fig, ax = plt.subplots(figsize=(15, 5))
    idxs = np.arange(len(target))
    for j, nucleotide in enumerate('ACGT'):
        loc_ddGs = []
        loc_ddG_error = []
        ticks = []
        for i, t in enumerate(target):
            seq = target[:i] + nucleotide + target[i+1:]
            if seq in seq_ddGs:
                loc_ddGs.append(seq_ddGs[seq]/1000.0)
                loc_ddG_error.append(seq_ddG_error[seq]/1000.0)
                ticks.append(idxs[i])
            else:
                log.debug(i, nucleotide, target[i])
        ticks = np.array(ticks)
        ax.bar(ticks - 0.25 + 0.5 * j / 4.0, loc_ddGs,
               width=0.125,
               yerr=loc_ddG_error,
               color=base_colors[nucleotide],
               error_kw=dict(ecolor='k', alpha=0.6),
               label=nucleotide)
    ax.xaxis.grid(False)
    ax.set_xlim((-0.5, len(target)-0.5))
    ax.set_xticks(range(len(target)))
    ax.set_xticklabels(target)

    ylim = ax.get_ylim()
    for i, nucleotide in enumerate(target):
        ax.fill_between([i-0.5, i+0.5], [ylim[0]]*2, [ylim[1]]*2,
                        color=base_colors[nucleotide],
                        alpha=0.07)
    ax.set_ylim(ylim)
    ax.set_title('Single Mismatch $\Delta \Delta G$\'s', fontsize=fs)
    ax.set_xlabel('{} Reference Sequence (Background Color)'.format(reference_sequence_name), fontsize=fs)
    ax.set_ylabel(r'$\Delta \Delta G \left(\frac{kJ}{mol} \right)$', fontsize=fs)
    ax.legend(loc='best')
    return fig


def plot_kd_list(seq_Kds, ham_seqs, hamming_distance):
    assert hamming_distance >= 0
    kd_list = np.array([seq_Kds[seq]/1000.0 for seq in ham_seqs if seq in seq_Kds])
    kd_list.sort()
    fig, ax = plt.subplots(figsize=(15, 7))
    ax.plot(range(len(kd_list)), kd_list)
    ax.set_yscale('log')
    ax.set_xlabel('HamDist=%d Seqs Sorted by $K_d$' % hamming_distance)
    ax.set_ylabel('$K_d$ (nM)')
    ax.set_title('HamDist=%d Sorted $K_d$\'s' % hamming_distance)
    return fig


def write_ddgs(ddGs, ddG_error, filename):
    # for fname, ddGs, ddG_error in [('target{}_close_seq_ddGs_and_errors.txt'.format(target_name), seq_ddGs, seq_ddG_error)]:
    with open(filename, 'w') as out:
        out.write('# Seq\tddG (J/mol)\tddG error (J/mol)\n')
        out.write('\n'.join(['%s\t%f\t%f' % (seq, ddGs[seq], ddG_error[seq]) for seq in sorted(ddGs.keys())]))


def write_kds(Kds, Kd_error, filename):
    with open(filename, 'w') as out:
        out.write('# Seq\tKd (pM)\tKd error (pM)\n')
        out.write('\n'.join(['%s\t%f\t%f' % (seq, Kds[seq], Kd_error[seq]) for seq in sorted(Kds.keys())]))


def plot_good_ham_reads(counts, bins, zoom_in, hamming_distance):
    if zoom_in:
        fig, ax = plt.subplots()
    else:
        fig, ax = plt.subplots(figsize=(15, 6))
    ax.hist(counts, bins, histtype='step')
    ax.set_title('Good Ham=%d Reads Found' % hamming_distance)
    ax.set_xlabel('Read Counts')
    ax.set_ylabel('Number of Seqs')
    if zoom_in:
        ax.set_xlim((0, 50))
    return fig


def delta_G(Kd):
    """Takes a Kd in pM and return delta_G in J/mol"""
    return constants.R * 333 * np.log(Kd / 10 ** 12)


def delta_delta_G(Kd, ref_delta_G):
    """Returns delta_delta_G in J/mol"""
    return delta_G(Kd) - ref_delta_G


def calculate_ddg(h5_filepaths, int_scores, close_reads, ref_delta_G, protein_channel, fobs_func):
    seq_Kds, seq_Kd_error, seq_ddGs, seq_ddG_error = {}, {}, {}, {}
    for seq, read_names in close_reads.items():
        if len(read_names) < 5:
            continue
        read_names = list(read_names)
        popt = curve_fit_Fobs_fixed_curve_given_read_names(int_scores, h5_filepaths, read_names, protein_channel, fobs_func)
        seq_Kds[seq] = popt[0]
        seq_ddGs[seq] = delta_delta_G(popt[0], ref_delta_G)

        bootstrap_Kds, bootstrap_ddGs = [], []
        for _ in range(50):
            resamp_read_names = np.random.choice(read_names, size=len(read_names), replace=True)
            try:
                popt = curve_fit_Fobs_fixed_curve_given_read_names(int_scores, h5_filepaths, resamp_read_names,
                                                                   protein_channel, fobs_func)
            except:
                log.error("%s, Read name length: %d, Resample read names length: %d" % (seq, len(read_names), len(resamp_read_names)))
                error.fail("Error calculating ddG")
            bootstrap_Kds.append(popt[0])
            bootstrap_ddGs.append(delta_delta_G(popt[0], ref_delta_G))
        seq_Kd_error[seq] = np.std(bootstrap_Kds)
        seq_ddG_error[seq] = np.std(bootstrap_ddGs)
    return seq_Kds, seq_Kd_error, seq_ddGs, seq_ddG_error


def plot_fluorescence_vs_concentration(intensities, kd, fmax, fmin, fobs_func, nM_concentrations):
    xx = np.logspace(1, 5, 200)
    yy = fobs_func(xx, kd)

    fig, ax = plt.subplots(figsize=(10, 7))
    for intensity in intensities:
        ax.plot(nM_concentrations, intensity, 'b', alpha=0.01)
    ax.plot(xx / 1000, yy, 'r', label='Fit Curve', linewidth=2.5)
    ax.set_xscale('log')
    ax.grid(False)

    ax.set_title('Fluorescence vs Concentration Curve')
    ax.set_xlabel('Concentration (nM)')
    ax.set_ylabel('Intensity')
    xlim, ylim = ax.get_xlim(), ax.get_ylim()

    inc = (ylim[1] - ylim[0]) / 3
    oom = int(np.log10(inc))
    inc -= inc % max(1, int(0.05 * 10 ** oom))
    ax.yaxis.set_major_locator(MultipleLocator(inc))

    ax.text(0.2,
            sum(ax.get_ylim()) / 2,
            '$K_d = %.2f$ nM\n$F_{max} = %.1f$\n$F_{min} = %.1f$' % (kd / 1000, fmax, fmin),
            fontsize=20,
            va='center')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(18)

    ax.set_axis_bgcolor('white')
    return fig


def curve_fit_Fobs_fixed_curve_given_read_names(int_scores, h5_filepaths, read_names, protein_channel, fobs_func):
    all_pM_concentrations = []
    all_intensities = []
    for h5_fpath in h5_filepaths:
        pM_conc = misc.parse_concentration(h5_fpath)
        if protein_channel not in int_scores.score_given_read_name_in_channel[h5_fpath]:
            continue
        score_dict = int_scores.score_given_read_name_in_channel[h5_fpath][protein_channel]
        for read_name in read_names:
            if read_name in score_dict:
                all_pM_concentrations.append(pM_conc)
                all_intensities.append(float(score_dict[read_name]))
    log.debug("all pM concentrations: {}".format(all_pM_concentrations))
    log.debug("all intensities: {}".format(all_intensities))
    return curve_fit(fobs_func, all_pM_concentrations, all_intensities)[0]


def fit_curve_given_read_names(int_scores, protein_channel, read_names, h5_filepaths):
    all_pM_concentrations = []
    all_intensities = []
    for h5_fpath in h5_filepaths:
        pM_conc = misc.parse_concentration(h5_fpath)
        if protein_channel not in int_scores.score_given_read_name_in_channel[h5_fpath]:
            continue
        score_dict = int_scores.score_given_read_name_in_channel[h5_fpath][protein_channel]
        for read_name in read_names:
            if read_name in score_dict:
                all_pM_concentrations.append(pM_conc)
                all_intensities.append(score_dict[read_name])
    optimization_result = minimize(make_Fobs_sq_error(h5_filepaths, all_pM_concentrations, all_intensities),
                                   (500, 1),
                                   bounds=((0, None), (0, None)))
    # we just return 'x', the solution array
    return optimization_result.x


def get_sequences_given_ref_and_hamming_distance(ref_seq, ham):
    seqs = []
    for idxs in itertools.combinations(range(len(ref_seq)), ham):
        mm_bases = ['ACGT'.replace(ref_seq[idx], '') for idx in idxs]
        for new_bases in itertools.product(*mm_bases):
            new_seq = ref_seq[:idxs[0]]
            for i, new_base in enumerate(new_bases[:-1]):
                new_seq += new_base + ref_seq[idxs[i]+1:idxs[i+1]]
            new_seq += new_bases[-1] + ref_seq[idxs[-1]+1:]
            seqs.append(new_seq)
    return seqs


def load_close_reads(read_names_by_seq_fpath, close_seqs, good_read_names):
    close_reads = {seq: set() for seq in close_seqs}
    for line in open(read_names_by_seq_fpath):
        words = line.strip().split()
        seq = words[0]
        read_names = words[1:]
        for close_seq in close_seqs:
            if close_seq in seq:
                close_reads[close_seq].update(rn for rn in read_names if rn in good_read_names)
                break
    return close_reads


def load_bad_read_names(read_names_by_seq_fpath, off_target, good_read_names):
    bad_read_names = set()
    for line in open(read_names_by_seq_fpath):
        words = line.strip().split()
        seq = words[0]
        read_names = words[1:]
        if off_target in seq:
            bad_read_names.update(rn for rn in read_names if rn in good_read_names)
    return bad_read_names


def get_fmin(h5_filepaths, protein_channel, int_scores, bad_read_names):
    for h5_fpath in h5_filepaths:
        if protein_channel not in int_scores.score_given_read_name_in_channel[h5_fpath]:
            continue
        score_dict = int_scores.score_given_read_name_in_channel[h5_fpath][protein_channel]
        intensities = []
        for read_name in bad_read_names:
            if read_name in score_dict:
                intensities.append(score_dict[read_name])
        if intensities:
            break
            # TODO: Move return statement here?
    return np.average(intensities)


def Fobs(x, Kd, Fmax, Fmin):
    return Fmax / (1.0 + (float(Kd)/x)) + Fmin


def make_Fobs_sq_error(concentrations, intensities, Fmin):
    def Fobs_sq_error(params):
        Kd, Fmax = params
        return sum((Fobs(conc, Kd, Fmax, Fmin) - obs_avg)**2 for conc, obs_avg in zip(concentrations, intensities))
    return Fobs_sq_error


def sort_h5_files(data_directory):
    h5_filepaths = glob.glob(os.path.join(data_directory, '*.h5'))
    i = 0
    while i < len(h5_filepaths):
        if 'phix' in h5_filepaths[i].lower() or 'chip' in h5_filepaths[i]:
            h5_filepaths.pop(i)
        else:
            i += 1
    h5_filepaths.sort(key=misc.parse_concentration)
    return h5_filepaths


def calculate_intensities(good_perfect_read_names, sample_size, int_scores, protein_channel, h5_filepaths):
    intensities = []
    for read_name in random.sample(good_perfect_read_names, sample_size):
        intensity = [int_scores.score_given_read_name_in_channel[h5_fpath][protein_channel][read_name]
                     if protein_channel in int_scores.score_given_read_name_in_channel[h5_fpath]
                         and read_name in int_scores.score_given_read_name_in_channel[h5_fpath][protein_channel]
                     else None
                     for h5_fpath in h5_filepaths]
        intensities.append(intensity)
    return intensities


def fob_fix(fmin, fmax, x, kd):
    return fmax / (1.0 + (float(kd) / x)) + fmin


def calculate_nM_concentrations(h5_filepaths):
    return [misc.parse_concentration(h5_fpath) / 1000.0 for h5_fpath in h5_filepaths]


def determine_protein_channels(image_directory, metadata):
    channels = projectinfo.load_channels(image_directory)
    alignment_channel = {[metadata['alignment_channel']]}
    return channels - alignment_channel


def main(metadata, image_directory, target_info):
    output_directory = functools.partial(os.path.join, 'figures')
    protein_channels = determine_protein_channels(image_directory, metadata)
    read_names_by_seq_fpath = os.path.join(metadata.read_directory, 'read_names_by_seq.txt')
    perfect_target_read_name_fpath = os.path.join(metadata.read_directory,
                                                  'perfect_target_{}_read_names.txt'.format(target_info.on_target_label))
    perfect_target_read_names = set(line.strip() for line in open(perfect_target_read_name_fpath))
    h5_filepaths = sort_h5_files(image_directory)
    results_dirs = [os.path.join(image_directory, os.path.splitext(os.path.basename(h5_fpath))[0])
                    for h5_fpath in h5_filepaths]

    log.debug('Loading data...')
    int_scores = IntensityScores(h5_filepaths)
    int_scores.get_LDA_scores(results_dirs, metadata['lda_weights'])
    log.debug('Normalizing data...')
    int_scores.normalize_scores()
    int_scores.plot_aligned_images('br', 'o*')
    int_scores.plot_normalization_constants()
    int_scores.print_reads_per_channel()
    good_num_ims_cutoff = len(h5_filepaths) - 1
    int_scores.build_good_read_names(good_num_ims_cutoff)
    good_read_names = int_scores.good_read_names
    good_perfect_read_names = perfect_target_read_names & good_read_names
    log.debug('Good Perfect Reads: %d' % len(good_perfect_read_names))

    single_ham_seqs = get_sequences_given_ref_and_hamming_distance(target_info.on_target_sequence, 1)
    double_ham_seqs = get_sequences_given_ref_and_hamming_distance(target_info.on_target_sequence, 2)
    close_seqs = [target_info.on_target_sequence] + single_ham_seqs + double_ham_seqs
    close_reads = load_close_reads(read_names_by_seq_fpath, close_seqs, good_read_names)
    single_counts = [len(close_reads[seq]) for seq in single_ham_seqs]
    double_counts = [len(close_reads[seq]) for seq in double_ham_seqs]
    int_scores.build_score_given_read_name_given_channel()

    bad_read_names = load_bad_read_names(read_names_by_seq_fpath, target_info.off_target_sequence, good_read_names)
    sample_size = min(2000, len(good_perfect_read_names))

    for protein_channel in protein_channels:
        Fmin = get_fmin(h5_filepaths, protein_channel, int_scores, bad_read_names)
        Kd, Fmax = fit_curve_given_read_names(int_scores,
                                              protein_channel,
                                              random.sample(good_perfect_read_names, sample_size),
                                              h5_filepaths)
        Fobs_fixed = functools.partial(fob_fix, Fmin, Fmax)
        nM_concentrations = calculate_nM_concentrations(h5_filepaths)
        ref_delta_G = delta_G(Kd)

        log.info("Good read names in {}: {}".format(protein_channel, len(good_read_names)))

        plot_good_ham_reads(single_counts, 50, False, 2).savefig("{}_good_ham_reads_50.png")
        plot_good_ham_reads(single_counts, 50, True, 2).savefig("{}_good_ham_reads_50_zoomed.png")
        plot_good_ham_reads(double_counts, 200, False, 2).savefig("{}_good_ham_reads_200.png")
        plot_good_ham_reads(double_counts, 200, True, 2).savefig("{}_good_ham_reads_200_zoomed.png")

        seq_Kds, seq_Kd_error, seq_ddGs, seq_ddG_error = calculate_ddg(h5_filepaths, int_scores, close_reads,
                                                                       ref_delta_G, protein_channel, Fobs_fixed)

        write_kds(seq_Kds, seq_Kd_error,
                  'target{}_{}_close_seq_Kds_and_errors.txt'.format(target_info.on_target_label, protein_channel))
        write_ddgs(seq_ddGs, seq_ddG_error,
                   'target{}_{}_close_seq_ddGs_and_errors.txt'.format(target_info.on_target_label, protein_channel))
        single_ham_seqs_fig = plot_kd_list(seq_Kds, single_ham_seqs, 1)
        single_ham_seqs_fig.savefig(output_directory('{}_single_ham_seqs.png'.format(protein_channel)))
        double_ham_seqs_fig = plot_kd_list(seq_Kds, double_ham_seqs, 2)
        double_ham_seqs_fig.savefig(output_directory('{}_double_ham_seqs.png'.format(protein_channel)))
        single_mismatch_ddgs_fig = plot_single_mismatch_ddgs(seq_ddGs, seq_ddG_error, target_info.on_target_sequence,
                                                             'Target %s' % target_info.on_target_label.upper(), fs=16)
        single_mismatch_ddgs_fig.savefig(output_directory('{}_single_mismatch_ddgs.png'.format(protein_channel)))

        intensities = calculate_intensities(good_perfect_read_names, sample_size, int_scores, protein_channel, h5_filepaths)
        fl_vs_conc_fig = plot_fluorescence_vs_concentration(intensities, Kd, Fmax, Fmin, Fobs_fixed, nM_concentrations)
        fl_vs_conc_fig.savefig(output_directory('{}_fl_vs_conc.png'.format(protein_channel)))
