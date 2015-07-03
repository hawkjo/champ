import sys
import math
import random
import time
from collections import Counter, defaultdict
from itertools import izip
from Bio import SeqIO
from general_sequence_tools import dna_rev_comp
from adapters_cython import simple_hamming_distance


def approx_eq(a, b, tol=0.001):
    return abs(a - b) < tol


def collapse_if_overlapped_pair(seq_list, min_overlap=10, frac_error=0.1):
    if len(seq_list) == 1:
        return seq_list[0]
    elif len(seq_list) == 2:
        R1 = seq_list[0]
        R2_rc = dna_rev_comp(seq_list[1])
        max_overlap = min(map(len, [R1, R2_rc]))
        for i in range(min_overlap, max_overlap):
            if simple_hamming_distance(R2_rc[:i], R1[-i:]) < frac_error * i:
                return R1 + R2_rc[i:]
            elif simple_hamming_distance(R1[:i], R2_rc[-i:]) < frac_error * i:
                return R2_rc + R1[i:]
    return tuple(seq_list)  # If none of the above


class GibbsMotifSampler:
    def __init__(self):
        self.bases = 'ACGT'

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
        tmp_seqs = defaultdict(list)
        for rec in SeqIO.parse(open(fpath), 'fastq'):
            tmp_seqs[rec.id].append(str(rec.seq))
        self.seqs = [collapse_if_overlapped_pair(seq_list) for seq_list in tmp_seqs.values()]
        self.num_seqs = len(self.seqs)

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
            return [s[i:i+self.k] for i in range(len(s) - self.k + 1)]
        else:
            # Other option is to introduce iterable of strings
            return [ss[i:i+self.k] for ss in s for i in range(len(ss) - self.k + 1)]

    def random_kmer(self, s):
        if isinstance(s, str):
            i = random.randrange(0, len(s) - self.k + 1)
            return s[i:i+self.k]
        else:
            # Using all_kmers incorporates the robustness found there
            return random.choice(self.all_kmers(s))

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

    def gibbs_sample(self):
        start_time = time.time()
        self.motifs = [self.random_kmer(seq) for seq in self.seqs]
        self.make_profile_counts()
        absolute_best_motifs = self.motifs[:]
        absolute_best_score = self.motifs_score()
        
        for i in range(20):
            print
            print i
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

    def print_motifs(self):
        print 'Motif score:', self.motifs_score()
        print '\n'.join(self.motifs)

if __name__ == '__main__':
    gms = GibbsMotifSampler()
    gms.parse_rosalind_file(sys.argv[1])
    gms.gibbs_sample()
    gms.print_motifs()
