from collections import defaultdict
import sys
from champ import adapters_cython
import scipy.misc
import matplotlib as mpl
import matplotlib.colors as mcolors
import numpy as np


def build_read_names_given_seq(target,
                               read_names_by_seq_fpath,
                               allowed_read_names_set,
                               is_interesting_seq,
                               max_ham,
                               verbose=True):
    interesting_reads = defaultdict(set)
    for i, line in enumerate(open(read_names_by_seq_fpath)):
        if verbose and i % 10000 == 0:
            sys.stdout.write('.')
            sys.stdout.flush()

        words = line.strip().split()
        seq = words[0]
        if is_interesting_seq(seq):
            read_names = set(words[1:]) & allowed_read_names_set
            interesting_reads[seq].update(read_names)
            last_start = len(seq) - len(target)
            if last_start < 0:
                continue
            min_ham_idx = min(range(0, last_start+1),
                              key=lambda i: adapters_cython.simple_hamming_distance(target, seq[i:i+len(target)]))
            min_ham = adapters_cython.simple_hamming_distance(target, seq[min_ham_idx:min_ham_idx+len(target)])
            if min_ham <= max_ham:
                min_ham_seq = seq[min_ham_idx:min_ham_idx+len(target)]
                interesting_reads[min_ham_seq].update(read_names)
    return interesting_reads


def plot_library_comp_by_hamming_distance(ax,
                                          target,
                                          max_ham,
                                          min_reads,
                                          interesting_reads,
                                          interesting_seqs):
    read_counts_given_ham_dist = defaultdict(list)
    for seq in interesting_seqs:
        if len(seq) != len(target):
            continue
        ham_dist = adapters_cython.simple_hamming_distance(target, seq)
        if ham_dist > max_ham:
            continue
        nreads = len(interesting_reads[seq])
        if nreads >= min_reads:
            read_counts_given_ham_dist[ham_dist].append(nreads)

    def nseqs_given_ham(ham_dist):
        return scipy.misc.comb(len(target), ham_dist) * 3 ** ham_dist

    def make_colormap(seq):
        """Return a LinearSegmentedColormap
        seq: a sequence of floats and RGB-tuples. The floats should be increasing
        and in the interval (0,1).
        """
        seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
        cdict = {'red': [], 'green': [], 'blue': []}
        for i, item in enumerate(seq):
            if isinstance(item, float):
                r1, g1, b1 = seq[i - 1]
                r2, g2, b2 = seq[i + 1]
                cdict['red'].append([item, r1, r2])
                cdict['green'].append([item, g1, g2])
                cdict['blue'].append([item, b1, b2])
        return mcolors.LinearSegmentedColormap('CustomMap', cdict)

    bar_w = 0.8
    viol_w = 0.5

    max_count = 0
    ham_dists, nseqs, log_fracs = [], [], []
    for ham_dist, good_count_list in sorted(read_counts_given_ham_dist.items()):
        frac_found = float(len(good_count_list)) / nseqs_given_ham(ham_dist)
        ham_dists.append(ham_dist)
        nseqs.append(len(good_count_list))
        log_fracs.append(np.log10(frac_found))
        max_count = max(good_count_list + [max_count])
    upper_xlim = 10 ** (int(np.log10(max_count)) + 2)
    text_x = 10 ** ((np.log10(upper_xlim) + np.log10(max_count)) / 2.0)

    min_log_frac = -5.0
    high_color = np.array([0.3, 0.3, 1])
    low_color = 0.8 * np.array([1, 1, 1])
    for ham_dist, nseqs, log_frac in zip(ham_dists, nseqs, log_fracs):
        if log_frac < min_log_frac:
            log_frac = min_log_frac
        color = low_color + (log_frac - min_log_frac) / (-min_log_frac) * (high_color - low_color)
        ax.barh(ham_dist - bar_w / 2.0, nseqs, height=bar_w, color=color, zorder=-1)
        ax.text(text_x, ham_dist, str(nseqs), ha='center', va='center')

    viol_color = 0.6 * np.array([1, 1, 1])
    for ham_dist, good_count_list in sorted(read_counts_given_ham_dist.items()):
        if len(good_count_list) > 1:
            viol_d = ax.violinplot(good_count_list, [ham_dist], showextrema=False, vert=False, widths=viol_w)
            viol_d['bodies'][0].set_color(viol_color)
            viol_d['bodies'][0].set_alpha(1)
            viol_d['bodies'][0].set_edgecolor('k')
        else:
            ax.plot([good_count_list[0]] * 2, [ham_dist - viol_w / 2.0, ham_dist + viol_w / 2.0], color='k', linewidth=1)

    ax.set_xscale('log')
    ax.set_xlim((0.7, upper_xlim))
    ax.plot([min_reads] * 2, ax.get_ylim(), ':k')

    ax.set_ylabel('Substitutions', fontsize=18)
    ax.set_yticks(range(max_ham + 1))
    ax.set_ylim((max_ham + 1, -1))
    ax.set_axis_bgcolor('white')
    ax.grid(False)
    for item in ax.get_xticklabels() + ax.get_yticklabels():
        item.set_fontsize(16)

    ax.set_xlabel('Unique Sequences (bars)\nClusters per Sequence (violins)', fontsize=18)

    cax, kw = mpl.colorbar.make_axes(ax)
    cmap = make_colormap([low_color, high_color])
    norm = mpl.colors.Normalize(vmin=min_log_frac, vmax=0)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cbar.set_label('Fraction of Sequences Recovered')
    ticks = range(int(min_log_frac), 1)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(['$\leq 10^{%d}$' % min_log_frac] + ['$10^{%d}$' % xx for xx in ticks[1:]])
    cbar.ax.tick_params(labelsize=14)
