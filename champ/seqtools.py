import sys
import itertools
import numpy as np
from collections import defaultdict
from champ.adapters_cython import simple_hamming_distance
import scipy.misc
import matplotlib as mpl
import matplotlib.colors as mcolors
from intensity import get_reasonable_process_count
from multiprocessing import Process
from multiprocessing.queues import SimpleQueue

bases = 'ACGT'


def mm_names(ref, seq):
    mms = []
    for i, (c1, c2) in enumerate(zip(ref, seq)):
        if c1 != c2:
            mms.append('{}{}{}'.format(c1, i, c2))
    return ','.join(mms)


def get_deletion_seqs(seq, ndel):
    """Returns set of all sequences with ndel deletions from given seq."""
    outset = set()
    for tup in itertools.combinations(range(len(seq)), r=ndel):
        newseq = seq[:tup[0]]
        for i, j in zip(tup, tup[1:]):
            newseq += seq[i + 1:j]
        newseq += seq[tup[-1] + 1:]
        assert len(newseq) == len(seq) - ndel, (tup, newseq)
        outset.add(newseq)
    return outset


def get_contiguous_insertion_seqs(seq, len_ins):
    """Returns set of all sequences with single insertions of length len_ins from given seq."""
    outset = set()
    all_insertions = [''.join(tup) for tup in itertools.product(bases, repeat=len_ins)]
    for i in range(1, len(seq) + 1):
        outset.update([seq[:i] + insertion + seq[i:] for insertion in all_insertions])
    assert all([len(outseq) == len(seq) + len_ins for outseq in outset])
    return outset


def get_insertion_seqs(seq, nins):
    """Returns set of all sequences with nins insertions from given seq."""
    outset = set()
    for tup in itertools.combinations(range(1, len(seq) + 1), r=nins):
        for ins_bases in itertools.product(bases, repeat=nins):
            assert len(ins_bases) == len(tup), (tup, ins_bases)
            newseq = seq[:tup[0]]
            for base_idx, (i, j) in enumerate(zip(tup, tup[1:])):
                newseq += ins_bases[base_idx] + seq[i:j]
            newseq += ins_bases[-1] + seq[tup[-1]:]
            assert len(newseq) == len(seq) + nins, (tup, newseq)
            outset.add(newseq)
    return outset


def get_mismatch_seqs(seq, num_mm):
    """Returns set of all sequences with num_mm mutations from given seq."""
    outset = set()
    for tup in itertools.combinations(range(len(seq)), r=num_mm):
        all_mm_bases = [bases.replace(seq[i], '') for i in tup]
        for mm_bases in itertools.product(*all_mm_bases):
            newseq = seq[:tup[0]]
            for i, c in enumerate(mm_bases[:-1]):
                newseq += c + seq[tup[i] + 1:tup[i + 1]]
            newseq += mm_bases[-1] + seq[tup[-1] + 1:]
            assert len(newseq) == len(seq), '{}\n{}'.format(seq, newseq)
            outset.add(newseq)
    return outset


complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def forward_complement(seq):
    return ''.join([complements[c] for c in seq])


def switch_end_to_complement(seq, num_bases):
    """Replaces num_bases bases of tail with complement"""
    if num_bases <= 0:
        return seq
    return seq[:-num_bases] + forward_complement(seq[-num_bases:])


def get_stretch_of_complement_seqs(seq, num_bases):
    """Returns all seqs with num_bases-long stretches of sequence replaced with complements"""
    outset = set()
    for i in range(len(seq) - num_bases + 1):
        outset.add(seq[:i] + forward_complement(seq[i:i + num_bases]) + seq[i + num_bases:])
    return outset


def get_randomized_stretch_seqs(seq, num_bases):
    """Returns all seqs with num_bases-long stretches of randomized nt."""
    outset = set()
    all_randomized = [''.join(tup) for tup in itertools.product(bases, repeat=num_bases)]
    for i in range(len(seq) - num_bases + 1):
        outset.update([seq[:i] + rand_seq + seq[i + num_bases:] for rand_seq in all_randomized])
    return outset


def get_randomized_pam_seqs(seq, num_pam_bases, num_randomized_bases, end='5p'):
    """Returns set of sequences with randomized pam and leading bases, at preferred end."""
    assert num_randomized_bases >= num_pam_bases
    all_randomized = (''.join(tup) for tup in itertools.product(bases, repeat=num_randomized_bases))
    if end == '5p':
        return set([rand_seq + seq[num_pam_bases:] for rand_seq in all_randomized])
    else:
        assert end == '3p', end
        return set([seq[:-num_pam_bases] + rand_seq for rand_seq in all_randomized])


def get_randomized_region_seqs(seq, start, end):
    """Returns set of sequences where seq[start:end] is randomized."""
    assert start < end, (start, end)
    all_randomized = (''.join(tup) for tup in itertools.product(bases, repeat=end - start))
    return set([seq[:start] + rand_seq + seq[end:] for rand_seq in all_randomized])


def get_mismatches_in_region(seq, start, end, num_mm):
    """Return all seqs with given number of mismatches in given region."""
    return set([seq[:start] + mm_seq + seq[end:] for mm_seq in get_mismatch_seqs(seq[start:end], num_mm)])


def get_complementary_bundle_sets(seq):
    """
    Return all sequences with combinations of stretches set to complementary sequence.

    For instance, take a sequence of length 13. Considering bundles of length 3 will
    produce the following bundles:

        ... ... ... ....

    Then forward-complimenting 2 bundles at a time will produce the following set of sequences:

        *** *** ... ....
        *** ... *** ....
        *** ... ... ****
        ... *** *** ....
        ... *** ... ****
        ... ... *** ****

    Note the last bundle includes left-over nucleotides.
    """
    outset = set()
    for bundle_len in range(3, 11, 2):  # only consider bundles up to length 10
        if bundle_len * 2 > len(seq):
            bundles = [(0, len(seq))]
        else:
            bundles = []
            for start in range(0, len(seq) - len(seq) % bundle_len - bundle_len, bundle_len):
                bundles.append((start, start + bundle_len))
            # Extend last bundle to end of sequence
            bundles.append((len(seq) - len(seq) % bundle_len - bundle_len, len(seq)))

        for num_bundles in range(2, 4):
            for selected_bundles in itertools.combinations(bundles, r=num_bundles):
                if len(seq) - sum(end - start for start, end in selected_bundles) <= len(seq) / 3.0:
                    # skip if there will be less than a third of the original sequence
                    continue
                newseq = seq[:]
                for start, end in selected_bundles:
                    newseq = newseq[:start] + forward_complement(newseq[start:end]) + newseq[end:]
                assert len(newseq) == len(seq)
                outset.add(newseq)
    return outset


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
            min_ham_idx = min(range(0, last_start + 1),
                              key=lambda i: simple_hamming_distance(target, seq[i:i + len(target)]))
            min_ham = simple_hamming_distance(target, seq[min_ham_idx:min_ham_idx + len(target)])
            if min_ham <= max_ham:
                min_ham_seq = seq[min_ham_idx:min_ham_idx + len(target)]
                interesting_reads[min_ham_seq].update(read_names)
    return interesting_reads


def _thread_build_interesting_sequences(read_name_sequences, interesting_sequences, results_queue):
    results = defaultdict(set)
    for rough_sequence, read_names in read_name_sequences:
        for interesting_sequence in interesting_sequences:
            if interesting_sequence in rough_sequence:
                results[interesting_sequence].update(read_names)
    results_queue.put(results)


def build_interesting_sequences(read_names_by_seq_filepath, interesting_sequences):
    process_count = get_reasonable_process_count()
    read_name_sequences = [[] for _ in range(process_count)]
    with open(read_names_by_seq_filepath) as f:
        for i, line in enumerate(f):
            words = line.strip().split()
            rough_sequence = words[0]
            read_names = set(words[1:])
            read_name_sequences[i % process_count].append((rough_sequence, read_names))

    results_queue = SimpleQueue()
    processes = []
    for rns in read_name_sequences:
        p = Process(target=_thread_build_interesting_sequences, args=(rns, interesting_sequences, results_queue))
        processes.append(p)
        p.start()

    interesting_read_names = defaultdict(set)
    for _ in read_name_sequences:
        results = results_queue.get()
        for interesting_sequence, read_names in results.items():
            interesting_read_names[interesting_sequence].update(read_names)
    for p in processes:
        p.join()
    return interesting_read_names


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
        ham_dist = simple_hamming_distance(target, seq)
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
