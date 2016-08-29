import glob
import itertools
import os
import random
import re
import sys
from collections import defaultdict
import h5py
import matplotlib.pyplot as plt
import numpy as np
from champ import misc, intensity, hdf5tools
from sklearn.neighbors import KernelDensity
import warnings
import yaml


def dot():
    sys.stdout.write(".")
    sys.stdout.flush()


def load_read_sequences(path):
    sequences = {}
    with open(path) as f:
        for line in f:
            data = line.split("\t")
            sequence = data[0]
            for read_name in data[1:]:
                sequences[read_name] = sequence
    return sequences


def load_target(name, target_file_path='/shared/targets.yml'):
    with open(target_file_path) as f:
        targets = yaml.load(f)
    return targets[name]


def load_read_name(path):
    with open(path) as f:
        return set(line.strip() for line in f)


def load_h5_filenames(directory, sort_by="concentration"):
    assert sort_by in ("concentration", "time", None)
    h5_paths = glob.glob(os.path.join(directory, '*.h5'))
    if sort_by == "concentration":
        h5_paths.sort(key=misc.parse_concentration)
    if sort_by == "time":
        raise NotImplementedError("Sorting by time is not yet supported.")
    return h5_paths


class Analysis(object):
    """
    Represents the various analyses that a scientist wants to perform on a dataset.
    It does not actually do any computation - instead, it just ensures that the right functions
    are run and that dependencies are satisfied.

    """
    def __init__(self, read_data_directory, image_directory):
        self._read_data_directory = read_data_directory
        self._image_directory = image_directory
        self._target_name = None
        self._target_sequence = None
        self._off_target_sequence = None
        self._lda_path = None
        self._read_name_paths = {}
        self._good_perfect_read_name_paths = set()
        self._analyses = set()
        self._plots = set()
        self.h5_paths = ()

    def add_read_names(self, name, filename):
        self._read_name_paths[name] = os.path.join(self._read_data_directory, filename)

    def analyze_hamming_distance(self):
        self._analyses.add("hamming_distance")

    def analyze_single_mismatch_penalties(self):
        self._analyses.add("single_mismatch")

    def analyze_kd(self, target_name, target_sequence, off_target_sequence, lda_path='/shared/bLDA_coef_nonneg.txt'):
        self._analyses.add("kd")
        self._target_name = target_name
        self._target_sequence = target_sequence
        self._off_target_sequence = off_target_sequence
        self._lda_path = lda_path

    def show_aligned_images(self):
        self._analyses.add("lda")
        self._plots.add("aligned_images")

    def show_normalization_constants(self):
        self._analyses.add("lda")
        self._plots.add("normalization_constants")

    @property
    def good_perfect_read_filenames(self):
        return {[os.path.join(self._read_data_directory, "perfect_target_{}_read_names.txt".format(self._target_name.lower())),
                 os.path.join(self._read_data_directory, "target_{}_read_names.txt".format(self._target_name.lower()))]}

    @property
    def lda_path(self):
        return self._lda_path

    @property
    def analyses(self):
        return self._analyses

    @property
    def results_directory(self):
        return os.path.join(self._image_directory, "results")

    @property
    def figure_directory(self):
        return os.path.join(self._image_directory, "figs")

    @property
    def results_directories(self):
        return [os.path.join(self.results_directory, os.path.splitext(os.path.basename(path))[0])
                for path in self.h5_paths]


def get_int_scores(sorted_h5_filepaths, results_directories, lda_path):
    int_scores = intensity.IntensityScores(sorted_h5_filepaths)
    int_scores.get_LDA_scores(results_directories, lda_path)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        int_scores.normalize_scores()
    return int_scores


def get_good_reads(int_scores, perfect_target_read_names, h5_filepaths):
    int_scores.print_reads_per_channel()
    good_num_ims_cutoff = len(h5_filepaths) - 1
    int_scores.build_good_read_names(good_num_ims_cutoff)
    good_read_names = int_scores.good_read_names
    good_perfect_read_names = perfect_target_read_names & good_read_names
    print('Good Perfect Reads:', len(good_perfect_read_names))
    return good_read_names, good_perfect_read_names


def hamming_distance(perfect_target_sequence, good_read_names, read_name_directory):
    bases = 'ACGT'

    def get_sequences_given_ref_and_hamming_distance(ref_seq, ham):
        seqs = []
        for idxs in itertools.combinations(range(len(ref_seq)), ham):
            mm_bases = [bases.replace(ref_seq[idx], '') for idx in idxs]
            for new_bases in itertools.product(*mm_bases):
                new_seq = ref_seq[:idxs[0]]
                for i, new_base in enumerate(new_bases[:-1]):
                    new_seq += new_base + ref_seq[idxs[i] + 1:idxs[i + 1]]
                new_seq += new_bases[-1] + ref_seq[idxs[-1] + 1:]
                seqs.append(new_seq)
        return seqs

    single_ham_seqs = get_sequences_given_ref_and_hamming_distance(perfect_target_sequence, 1)
    double_ham_seqs = get_sequences_given_ref_and_hamming_distance(perfect_target_sequence, 2)
    close_seqs = [perfect_target_sequence] + single_ham_seqs + double_ham_seqs

    close_reads = {seq: set() for seq in close_seqs}
    read_names_by_seq_fpath = os.path.join(read_name_directory, 'read_names_by_seq.txt')
    for line in open(read_names_by_seq_fpath):
        words = line.strip().split()
        seq = words[0]
        read_names = words[1:]
        for close_seq in close_seqs:
            if close_seq in seq:
                close_reads[close_seq].update(rn for rn in read_names if rn in good_read_names)
                break

    single_counts = [len(close_reads[seq]) for seq in single_ham_seqs]
    double_counts = [len(close_reads[seq]) for seq in double_ham_seqs]

    fig, ax = plt.subplots(figsize=(15, 6))
    ax.hist(single_counts, 50, histtype='step')
    ax.set_title('Good Ham=1 Reads Found')
    ax.set_xlabel('Read Counts')
    ax.set_ylabel('Number of Seqs')

    fig, ax = plt.subplots()
    ax.hist(single_counts, 50, histtype='step')
    ax.set_title('Good Ham=1 Reads Found')
    ax.set_xlabel('Read Counts')
    ax.set_ylabel('Number of Seqs')
    ax.set_xlim((0, 50))
    plt.show()
    plt.close()

    fig, ax = plt.subplots(figsize=(15, 6))
    ax.hist(double_counts, 200, histtype='step')
    ax.set_title('Good Ham=2 Reads Found')
    ax.set_xlabel('Read Counts')
    ax.set_ylabel('Number of Seqs')

    fig, ax = plt.subplots()
    ax.hist(double_counts, 200, histtype='step')
    ax.set_title('Good Ham=2 Reads Found')
    ax.set_xlabel('Read Counts')
    ax.set_ylabel('Number of Seqs')
    ax.set_xlim((0, 50))
    plt.show()
    plt.close()


def plot_aligned_images(int_scores):
    for _ in int_scores.plot_aligned_images('br', 'o*'):
        plt.show()
    plt.close()


def plot_normalization_constants(int_scores):
    for _ in int_scores.plot_normalization_constants():
        plt.show()
    plt.close()


def validate(analysis):
    # checks if all external files needed for an analysis exist
    pass


def run(analysis):
    pass



# ######################################################################3
# Example!
target_sequence = load_target("D")
off_target_sequence = load_target("B")
analysis = Analysis('.', '/shared/SA16105/read_names')
analysis.analyze_hamming_distance()
analysis.analyze_single_mismatch_penalties()
analysis.analyze_kd("D", target_sequence, off_target_sequence)
validate(analysis)
run(analysis)
