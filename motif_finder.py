import sys
import math
import random
import time
from copy import deepcopy
from collections import Counter, defaultdict
from itertools import izip
from Bio import SeqIO
from general_sequence_tools import dna_rev_comp
import matplotlib.pyplot as plt
import numpy as np
from fastqtools import collapse_if_overlapped_pair


def approx_eq(a, b, tol=0.001):
    return abs(a - b) < tol


class MotifFinder:
    def __init__(self):
        self.bases = 'ACGT'
        self.ignore_kmers = set()

    def parse_rosalind_file(self, fpath):
        with open(fpath) as f:
            var = f.readline().strip().split()
            self.k = int(var[0])  # length of motif
            self.num_seqs = int(var[1])  # num sequences
            self.num_iters = int(var[2])
            self.seqs = []
            line = f.readline().strip()
            while line:
                self.seqs.append(line.upper())
                line = f.readline().strip()
            assert self.num_seqs == len(self.seqs)

    def parse_fastq_file(self, fpath):
        self.input_biopython_records(SeqIO.parse(open(fpath), 'fastq'))

    def input_biopython_records(self, records):
        tmp_seqs = defaultdict(list)
        for rec in records:
            tmp_seqs[rec.id].append(str(rec.seq))
        self.set_seqs(tmp_seqs.values())

    def set_seqs(self, seqs):
        self.seqs = [collapse_if_overlapped_pair(seq_list) for seq_list in seqs]
        self.num_seqs = len(self.seqs)

    def set_ignore_kmers(self, ignore_strings):
        self.ignore_kmers = set([kmer for s in ignore_strings for kmer in self.all_kmers(s)])

    def set_params(self, **kwargs):
        for k, v in kwargs.items():
            assert k in ['k', 'num_iters'], k
            setattr(self, k, int(v))

    def make_profile_counts(self):
        for motif in self.motifs:
            assert len(motif) == len(self.motifs[0]), 'Unequal motifs'
        self.profile_counts = {}
        for c in self.bases:
            self.profile_counts[c] = [1 + sum(1 for motif in self.motifs if motif[i] == c)
                                      for i in range(self.k)]  # Laplace smoothing used
        for i in range(self.k):
            assert sum(self.profile_counts[c][i] for c in self.bases) == self.num_seqs + len(self.bases)
    
    def make_profile(self):
        # The denominator isn't as simple as it could be because sometimes we will be omitting one
        # motif from the bunch.
        denom = float(sum(self.profile_counts[c][0] for c in self.bases))
        self.profile = {c: [self.profile_counts[c][i] / denom for i in range(self.k)]
                        for c in self.bases}
        for i in range(len(self.motifs[0])):
            assert approx_eq(sum(self.profile[c][i] for c in self.bases), 1), 'Non-uniform probability'

    def motifs_score(self):
        score = 0
        for i in range(self.k):
            consensus_base = max(self.bases, key=lambda c: self.profile_counts[c][i])
            score += sum(self.profile_counts[c][i] for c in self.bases if c != consensus_base)
        return score

    def pattern_prob(self, pattern):
        return math.exp(sum(math.log(self.profile[c][i]) for i, c in enumerate(pattern)))

    def all_kmers(self, s):
        """All substrings of a given string of given length."""
        if isinstance(s, str):
            return [ss[i:i+self.k]
                    for ss in [s, dna_rev_comp(s)]
                    for i in range(len(ss) - self.k + 1)
                    if ss[i:i+self.k] not in self.ignore_kmers]
        else:
            # Other option is to introduce iterable of strings
            it = deepcopy(s)
            return [ss[i:i+self.k]
                    for s in it
                    for ss in [s, dna_rev_comp(s)]
                    for i in range(len(ss) - self.k + 1)
                    if ss[i:i+self.k] not in self.ignore_kmers]

    def random_kmer(self, s):
        return random.choice(self.all_kmers(s))

    def print_motifs(self):
        print 'Motif score:', self.motifs_score()
        print '\n'.join(self.motifs)


class GibbsMotifSampler(MotifFinder):
    def profile_random_motif(self, i):
        substrings = self.all_kmers(self.seqs[i])
        probs = [self.pattern_prob(motif) for motif in substrings]
        sum_probs = sum(probs)
        r = random.uniform(0, sum_probs)
        upto = 0
        for s, p in izip(substrings, probs):
            if upto + p > r:
                return s
            upto += p
        assert False, 'Shouldn\'t get here'

    def get_new_motif(self):
        idx = random.randrange(0, self.num_seqs)

        old_motif = self.motifs[idx]
        for i, c in enumerate(old_motif):
            self.profile_counts[c][i] -= 1
            assert self.profile_counts[c][i] != 0, 'Smoothed profile count down to 0.'
        self.make_profile()

        new_motif = self.profile_random_motif(idx)
        for i, c in enumerate(new_motif):
            self.profile_counts[c][i] += 1
        self.make_profile()

        self.motifs[idx] = new_motif

    def gibbs_sample(self, num_reps=10):
        start_time = time.time()
        self.motifs = [self.random_kmer(seq) for seq in self.seqs]
        self.make_profile_counts()
        absolute_best_motifs = self.motifs[:]
        absolute_best_score = self.motifs_score()
        
        for i in range(num_reps):
            self.motifs = [self.random_kmer(seq) for seq in self.seqs]
            best_motifs = self.motifs[:]
            best_score = self.motifs_score()
            self.make_profile_counts()
            for j in range(self.num_iters):
                if j % 1000 == 0:
                    sys.stdout.write('.')
                    sys.stdout.flush()
                self.get_new_motif()
                score = self.motifs_score()
                if score < best_score:
                    best_motifs = self.motifs[:]
                    best_score = score
            if best_score < absolute_best_score:
                absolute_best_motifs = best_motifs[:]
                absolute_best_score = best_score

        self.motifs = absolute_best_motifs[:]
        self.make_profile_counts()
        self.make_profile()

        print
        print 'Gibbs Sampler Time:', time.time() - start_time


class CommonKmerFinder(MotifFinder):
    def count_and_sort_kmers(self):
        self.kmer_counter = Counter()
        for seq in self.seqs:
            self.kmer_counter.update(set(self.all_kmers(seq)))  # set to count unique seqs
        self.kmers_and_counts = self.kmer_counter.items()
        self.kmers_and_counts.sort(key=lambda tup: tup[1], reverse=True)

    def print_n_most_common_kmers(self, n=30):
        print '\n'.join('%s:  %d' % (kmer, count) for kmer, count in self.kmers_and_counts[:n])

    def plot_n_most_common_kmers(self, n=20, ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 7))
        kmers = [kmer for kmer, count in self.kmers_and_counts[:n]]
        counts = [count for kmer, count in self.kmers_and_counts[:n]]
        pos = range(n)[::-1]
        ax.barh(pos, counts, align='center')
        plt.yticks(pos, kmers)
        ax.set_ylim((min(pos)-1, max(pos)+1))
        ax.set_xlabel('Counts')
        ax.set_title('Most Common %d-mers' % self.k)
        return ax


class EnrichedKmerFinder:
    def __init__(self, k, active_fpath, inactive_fpath):
        self.k = k

        self.active_ckm = CommonKmerFinder()
        self.active_ckm.parse_fastq_file(active_fpath)
        self.active_ckm.set_params(k=k)
        self.active_ckm.count_and_sort_kmers()

        self.inactive_ckm = CommonKmerFinder()
        self.inactive_ckm.parse_fastq_file(inactive_fpath)
        self.inactive_ckm.set_params(k=k)
        self.inactive_ckm.count_and_sort_kmers()

        self.kmers_not_in_inactive = set()
        self.fold_enrichment = {}
        for kmer, count in self.active_ckm.kmer_counter.items():
            if kmer in self.inactive_ckm.kmer_counter:
                self.fold_enrichment[kmer] = float(count) / self.inactive_ckm.kmer_counter[kmer]
            else:
                self.kmers_not_in_inactive.add(kmer)

    def print_n_most_common_kmers(self, n=20):
        print 'Top Fold Enrichments:'
        print '\n'.join('%s: %.3f' % (kmer, fe) for kmer, fe in
                        list(sorted(self.fold_enrichment.items(), key=lambda tup: -tup[1]))[:n])
        if self.kmers_not_in_inactive:
            print
            print 'Common kmers not in inactive sequences and counts:'
            print '\n'.join('%s: %d' % (kmer, count) for kmer, count in self.active_ckm.kmers_and_counts
                            if kmer in self.kmers_not_in_inactive and count > 10)

    def plot_n_most_common_kmers(self, n=20, ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 7))
        kmers = []; fes = []
        for kmer, fe in list(sorted(self.fold_enrichment.items(), key=lambda tup: -tup[1]))[:n]:
            kmers.append(kmer)
            fes.append(fe)
        pos = range(len(fes))[::-1]
        ax.barh(pos, fes, align='center')
        plt.yticks(pos, kmers)
        ax.set_ylim((min(pos)-1, max(pos)+1))
        ax.set_xlabel('Fold Enrichment')
        ax.set_title('Most Enriched %d-mers' % self.k)
        return ax


if __name__ == '__main__':
    gms = GibbsMotifSampler()
    gms.parse_rosalind_file(sys.argv[1])
    gms.gibbs_sample()
    gms.print_motifs()
