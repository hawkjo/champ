from collections import defaultdict
import numpy as np
from Bio.Seq import Seq


class TargetSequence(object):
    def __init__(self, sequence, pam_side=None, pam_length=None):
        """
        :param sequence: 
        :param pam_side: should be the integer 5 or 3 
        :param pam_length: 
        """
        self._sequence = sequence
        self._pam_side = int(pam_side)
        self._pam_length = pam_length
        self._human_indexes = {}

        if pam_length is not None:
            if self._pam_side == 3:
                self._guide_sequence = sequence[:-pam_length]
                for i in range(pam_length):
                    self._human_indexes[len(sequence) - pam_length + i] = -i - 1
                for i in range(len(sequence) - pam_length):
                    self._human_indexes[i] = len(sequence) - i - pam_length
            elif self._pam_side == 5:
                self._guide_sequence = sequence[pam_length:]
                for i in range(pam_length):
                    self._human_indexes[i] = i - pam_length
                for i in range(len(sequence) - pam_length):
                    self._human_indexes[i + pam_length] = i + 1
            else:
                raise ValueError("Invalid PAM side. Only use 3 or 5")
        else:
            for i in range(len(sequence)):
                self._human_indexes[i] = len(sequence) - i if pam_side == 3 else i + 1

    @property
    def pam(self):
        if self._pam_side == 3:
            return self._sequence[-self._pam_length:]
        return self._sequence[:self._pam_length]

    @property
    def pam_side(self):
        return self._pam_side

    @property
    def guide(self):
        return TargetSequence(self._guide_sequence, pam_side=self._pam_side)

    @property
    def sequence(self):
        return self._sequence

    @property
    def human_readable_indexes(self):
        return [(base, self._human_indexes[i]) for i, base in enumerate(self._sequence)]

    @property
    def single_deletions(self):
        for i in range(len(self._sequence)):
            yield i, self._sequence[:i] + self._sequence[i + 1:]

    @property
    def double_deletions(self):
        for i in range(len(self._sequence)):
            for j in range(i):
                seq = self._sequence[:j] + self._sequence[j + 1:i] + self._sequence[i + 1:]
                yield i, j, seq

    @property
    def single_mismatches(self):
        # produces all mismatch bases. Includes the correct base since we want to graph it
        for i in range(len(self._sequence)):
            for j, mismatch_base in enumerate('ACGT'):
                seq = self._sequence[:i] + mismatch_base + self._sequence[i + 1:]
                yield i, j, mismatch_base, seq

    @property
    def double_mismatches(self):
        bases = 'ACGT'
        # Double mismatch sequences
        for i in range(len(self._sequence)):
            for j in range(i):
                for mismatch_base_j in bases.replace(self._sequence[j], ''):
                    for mismatch_base_i in bases.replace(self._sequence[i], ''):
                        sequence = self._sequence[:j] + mismatch_base_j + self._sequence[j + 1:i] + mismatch_base_i + self._sequence[i + 1:]
                        yield i, j, mismatch_base_i, mismatch_base_j, sequence

        # Add in single mismatch data for comparison
        for i in range(len(self._sequence)):
            for mismatch_base in bases.replace(self._sequence[i], ''):
                sequence = self._sequence[:i] + mismatch_base + self._sequence[i + 1:]
                yield i, i, mismatch_base, mismatch_base, sequence

    @property
    def single_insertions(self):
        bases = 'ACGT'
        for j, insertion_base in enumerate(bases):
            for i, nt in enumerate(self._sequence):
                sequence = self._sequence[:i] + insertion_base + self._sequence[i:]
                yield i, j, insertion_base, sequence

    @property
    def double_insertions(self):
        bases = 'ACGT'
        for i in range(len(self._sequence)):
            for j in range(i):
                for insertion_base_1 in bases:
                    for insertion_base_2 in bases:
                        yield i, j, insertion_base_1, insertion_base_2, self._sequence[:j] + insertion_base_1 + self._sequence[j:i] + insertion_base_2 + self._sequence[i:]

        # single insertions for the diagonal
        for i in range(len(self._sequence)):
            for insertion_base in bases:
                yield i, i, insertion_base, insertion_base, self._sequence[:i] + insertion_base + self._sequence[i:]

    @property
    def complement_stretches(self):
        for stop in range(len(self._sequence)):
            for start in range(stop):
                yield start, stop, self._sequence[:start] + str(Seq(self._sequence[start:stop + 1]).complement()) + self._sequence[stop + 1:]


class TwoDMatrix(object):
    """
    Base class for storing 2D penalty data (where the axes are sequences of some kind)
    This was intended for mismatches and insertions, but now that I look at it, it might also work for deletions

    """

    def __init__(self, sequence, slots_per_position, bases='ACGT'):
        self._bases = bases
        self._sequence = sequence
        self._slots = slots_per_position
        self._values = defaultdict(dict)

    @property
    def _dimension(self):
        return self._slots * len(self._sequence)

    def to_matrix(self, side='lower', include_diagonal_values=True, flip_sequence=False):
        assert side in ('lower', 'upper')
        data = np.zeros((self._dimension, self._dimension))
        data[:] = np.nan
        for row, column_data in self._values.items():
            for column, value in column_data.items():
                if not include_diagonal_values and row == column:
                    continue
                if flip_sequence:
                    row = self._dimension - row - 1 + column
                    column = self._dimension - column - 1
                if side == 'lower':
                    data[row, column] = value
                elif side == 'upper':
                    data[column, row] = value
        return data


class MismatchMatrix(TwoDMatrix):
    def __init__(self, sequence, bases='ACGT'):
        super(MismatchMatrix, self).__init__(sequence, 3, bases)

    def set_value(self, position1, position2, base1, base2, value):
        index1 = self._bases.replace(self._sequence[position1], '').index(base1)
        index2 = self._bases.replace(self._sequence[position2], '').index(base2)
        r, c = position1 * self._slots + index1, position2 * self._slots + index2
        self._values[r][c] = value


class InsertionMatrix(TwoDMatrix):
    def __init__(self, sequence, bases='ACGT'):
        super(InsertionMatrix, self).__init__(sequence, 4, bases)

    def set_value(self, position1, position2, base1, base2, value):
        r, c = position1 * self._slots + self._bases.index(base1), position2 * self._slots + self._bases.index(base2)
        self._values[r][c] = value

    @property
    def data(self):
        return self._values


class DeletionMatrix(TwoDMatrix):
    def __init__(self, sequence):
        super(DeletionMatrix, self).__init__(sequence, 1, 'ACGT')

    def set_value(self, position1, position2, value):
        self._values[position1][position2] = value
