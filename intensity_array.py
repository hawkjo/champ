import numpy as np
import os
import re
import misc

bases = 'ACGT'
bases_set = set(bases)

class IntensityArray(object):
    def __init__(self):
        pass

    def parse_intensities_file(self, fpath):
        """
        Parses input file to read:
            course_trait_name
            course_trait_list
            h5_fpaths
            channel
            various attributes
            seqs
            read_names lol
            itensity_lolol 
                - list of lists of lists with intensity by seq by concentration by read.
                    Value of None for missing data

        build_derived_objects also called.
        """
        with open(fpath) as f:
            line = next(f)
            assert line.startswith('# Defining Course Trait:'), line
            self.course_trait_name = line.strip().split(': ')[1]
            self.course_trait_list = map(float, next(f).strip().split('\t'))
            line = next(f)
            assert line == '# HDF5 Files\n', line
            self.h5_fpaths = [next(f).strip() for i in range(len(self.course_trait_list))]
            for fpath in self.h5_fpaths:
                assert os.path.exists(fpath), fpath
            line = next(f)
            assert line.startswith('# Channel:'), line
            self.channel = line.strip().split(': ')[1]
            self.attr_names = []
            while True:
                line = next(f)
                if not line.startswith('#'):
                    break
                m = re.match('^# (.*): (.*)$', line)
                assert m, line
                setattr(self, m.group(1), m.group(2))
                self.attr_names.append(m.group(1))
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
                for _ in xrange(len(self.course_trait_list)):
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
        self.build_derived_objects()

    def build_derived_objects(self):
        """
        Sets derived traits, including:
            course_len
            intensity_loloarr
                - list of lists of arrays. Same as lolol, but last is np arrays, Nones removed
            idx_given_seq
            read_names_given_seq
            intensity_lol_given_seq
            intensity_loarr_given_seq
            nseqs
        """
        self.course_len = len(self.course_trait_list)
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

    def subIA(self, seqs=None, course_traits=None, max_clust=None):
        """
        Create an IntensityArray which is a subset of self with reduced seqs, course_traits, and/or
        max_clust.
        """
        assert seqs or course_traits or max_clust
        IA = IntensityArray()
        IA.course_trait_name = self.course_trait_name

        # Optionally reduce course_trait_list
        if course_traits:
            assert set(course_traits) <= set(self.course_trait_list), course_traits
            assert isinstance(course_traits, list)
            IA.course_trait_list = course_traits
            trait_idxs = [self.course_trait_list.index(trait) for trait in IA.course_trait_list]
            IA.h5_fpaths = [self.h5_fpaths[idx] for idx in trait_idxs]
        else:
            IA.course_trait_list = self.course_trait_list
            IA.h5_fpaths = self.h5_fpaths
            trait_idxs = range(len(IA.course_trait_list))

        # Optionally reduce seqs
        if seqs:
            IA.seqs = list(seqs)
        else:
            IA.seqs = self.seqs

        # Build intensity_lolol given reduced parameters
        IA.read_names = [self.read_names_given_seq[seq][:max_clust] for seq in IA.seqs]
        IA.intensity_lolol = []
        for seq in IA.seqs:
            old_lol = self.intensity_lol_given_seq[seq]
            IA.intensity_lolol.append([old_lol[idx][:max_clust] for idx in trait_idxs])

        # Copy attributes
        IA.channel = self.channel
        IA.attr_names = self.attr_names
        for attr_name in self.attr_names:
            setattr(IA, attr_name, getattr(self, attr_name))

        # Set derived 
        IA.build_derived_objects()
        return IA

    def medians_given_seq(self, seq):
        return map(np.median, self.intensity_loarr_given_seq[seq])

    def modes_given_seq(self, seq):
        return map(misc.get_mode, self.intensity_loarr_given_seq[seq])

    def stdevs_given_seq(self, seq):
        return map(np.std, self.intensity_loarr_given_seq[seq])

    def all_trait_and_inten_vals_given_seq(self, seq, max_clust=None):
        """
        Returns all concentration/intensities pairs in read_names, return as two lists.

        Returns:
            :list: all_concentrations
            :list: all_intensities
        """
        all_trait_vals, all_intensities = [], []
        for tval, inten_arr in zip(self.course_trait_list, self.intensity_loarr_given_seq[seq]):
            tmp_inten = list(inten_arr[:max_clust])
            all_trait_vals.extend([tval]*len(tmp_inten))
            all_intensities.extend(tmp_inten)
        return all_trait_vals, all_intensities

    def all_normalized_trait_and_inten_vals_given_seq(self, seq, Imin, Imax, max_clust=None):
        """
        Returns all concentration/intensities pairs in read_names, return as two lists, with the
        intensities adjusted by Imin and Imax to run typically between 0 and 1.

        Returns:
            :list: all_concentrations
            :list: all_intensities
        """
        def list_if_scalar(x):
            try:
                float(x)
                return [x]*len(self.course_trait_list)
            except:
                return x
        Imin = list_if_scalar(Imin)
        Imax = list_if_scalar(Imax)
        assert len(Imin) == len(Imax) == len(self.course_trait_list), (Imin, Imax)
        all_trait_vals, all_intensities = [], []
        for tval, imn, imx, inten_arr in zip(self.course_trait_list,
                                               Imin,
                                               Imax,
                                               self.intensity_loarr_given_seq[seq]):
            tmp_inten = list((inten_arr[:max_clust] - imn)/(imx - imn))
            all_trait_vals.extend([tval]*len(tmp_inten))
            all_intensities.extend(tmp_inten)
        return all_trait_vals, all_intensities
