import numpy as np
import os
import re

bases = 'ACGT'
bases_set = set(bases)

class IntensityArray(object):
    def __init__(self):
        pass

    def parse_intensities_file(self, fpath):
        """
        Parses input file.
        
        Results:
            course_trait_name
            course_trait_list
            course_len
            h5_fpaths
            channel
            various attributes
            seqs
            idx_given_seq
            itensity_lolol 
                - list of lists of lists with intensity by seq by concentration by read.
                    Value of None for missing data
            intensity_loloarr
                - list of lists of arrays. Same as lolol, but last is np arrays, Nones removed
        """
        with open(fpath) as f:
            line = next(f)
            assert line.startswith('# Defining Course Trait:'), line
            self.course_trait_name = line.strip().split(': ')[1]
            self.course_trait_list = map(float, next(f).strip().split('\t'))
            self.course_len = len(self.course_trait_list)
            line = next(f)
            assert line == '# HDF5 Files\n', line
            self.h5_fpaths = [next(f).strip() for i in range(self.course_len)]
            for fpath in self.h5_fpaths:
                assert os.path.exists(fpath), fpath
            line = next(f)
            assert line.startswith('# Channel:'), line
            self.channel = line.strip().split(': ')[1]
            while True:
                line = next(f)
                if not line.startswith('#'):
                    break
                m = re.match('^# (.*): (.*)$', line)
                assert m, line
                setattr(self, m.group(1), m.group(2))
            self.seqs = []
            self.read_names = []
            self.intensity_lolol = []
            while True:
                seq = line.strip()
                assert set(seq) <= bases_set, seq
                self.seqs.append(seq)
                line = next(f)
                self.read_names.append(line.strip().split('\t'))
                self.intensity_lolol.append([])
                for _ in xrange(self.course_len):
                    line = next(f)
                    self.intensity_lolol[-1].append(
                        [float(v) if v != '-' else None for v in line.strip().split('\t')]
                    )
                val_lens = map(len, self.intensity_lolol[-1])
                assert all(v == len(self.read_names[-1]) for v in val_lens), '\n'.join([seq] + self.intensity_lolol[-1])
                try:
                    line = next(f)
                except StopIteration:
                    break
            self.intensity_loloarr = []
            for seq_inten_list in self.intensity_lolol:
                self.intensity_loloarr.append([])
                for inten_list in seq_inten_list:
                    self.intensity_loloarr[-1].append(
                        np.array([v for v in inten_list if v is not None])
                    )
            self.idx_given_seq = {seq: i for i, seq in enumerate(self.seqs)}
            self.read_names_given_seq = {seq: self.read_names[i] for i, seq in enumerate(self.seqs)}
            self.intensity_lol_given_seq = {seq: self.intensity_lolol[i] for i, seq in enumerate(self.seqs)}
            self.intensity_loarr_given_seq = {seq: self.intensity_loloarr[i] for i, seq in enumerate(self.seqs)}
            self.nseqs = len(self.seqs)
