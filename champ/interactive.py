import glob
import re
from collections import defaultdict

import numpy as np
import yaml
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
        for downstream_index in range(len(self._sequence)):
            for upstream_index in range(downstream_index):
                seq = self._sequence[:upstream_index] + self._sequence[upstream_index + 1:downstream_index] + self._sequence[downstream_index + 1:]
                yield upstream_index, downstream_index, seq
            seq = self._sequence[:downstream_index] + self._sequence[downstream_index + 1:]
            yield downstream_index, downstream_index, seq

    @property
    def single_mismatches(self):
        # produces all mismatch bases. Includes the correct base since we want to graph it
        for position in range(len(self._sequence)):
            for base_index, mismatch_base in enumerate('ACGT'):
                seq = self._sequence[:position] + mismatch_base + self._sequence[position + 1:]
                yield position, base_index, mismatch_base, seq

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
        for downstream_index in range(len(self._sequence)):
            for upstream_index in range(downstream_index):
                for upstream_insertion_base in bases:
                    for downstream_insertion_base in bases:
                        yield downstream_index, upstream_index, upstream_insertion_base, downstream_insertion_base, self._sequence[:upstream_index] + upstream_insertion_base + self._sequence[upstream_index:downstream_index] + downstream_insertion_base + self._sequence[downstream_index:]

        # single insertions for the diagonal
        for i in range(len(self._sequence)):
            for insertion_base in bases:
                yield i, i, insertion_base, insertion_base, self._sequence[:i] + insertion_base + self._sequence[i:]

    @property
    def complement_stretches(self):
        for stop in range(len(self._sequence)):
            for start in range(stop):
                yield start, stop, self._sequence[:start] + str(Seq(self._sequence[start:stop + 1]).complement()) + self._sequence[stop + 1:]
        for position in range(len(self._sequence)):
            yield position, position, self._sequence[:position] + str(Seq(self._sequence[position]).complement()) + self._sequence[position + 1:]


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

    def to_matrix(self, side='lower', include_diagonal_values=True, flip_sequence=False, normalize_by=None):
        assert side in ('lower', 'upper')
        data = np.zeros((self._dimension, self._dimension))
        data[:] = np.nan
        for row, column_data in self._values.items():
            for column, values in column_data.items():
                if not include_diagonal_values and row == column:
                    continue
                if flip_sequence:
                    c = self._dimension - column - 1
                    r = self._dimension - row - 1
                else:
                    c = column
                    r = row
                if type(values) is list:
                    # we have multiple values for a single position, so we need to take the average in order to make
                    # a meaningful plot
                    clean_values = tuple(v for v in values if v is not None)
                    value = np.mean(clean_values) if clean_values else None
                else:
                    # values is just a single float, so alias it for the next few lines
                    value = values
                if normalize_by is not None and value is not None:
                    value /= normalize_by
                try:
                    if (side == 'lower' and not flip_sequence) or (side == 'upper' and flip_sequence):
                        data[r, c] = value
                    elif (side == 'upper' and not flip_sequence) or (side == 'lower' and flip_sequence):
                        data[c, r] = value
                except IndexError:
                    # In case the sequences aren't the same length. This is not a good solution
                    continue
        return data

    def _safe_append(self, r, c, value):
        current = self._values[r].get(c)
        if not current:
            self._values[r][c] = [value]
        else:
            self._values[r][c].append(value)


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
        r = position1 * self._slots + self._bases.index(base1)
        c = position2 * self._slots + self._bases.index(base2)
        self._values[r][c] = value


class SinglePositionMatrix(TwoDMatrix):
    """ Used for deletions, and for comparing two incompatible sequences by position alone. """
    def __init__(self, sequence):
        super(SinglePositionMatrix, self).__init__(sequence, 1, 'ACGT')

    def set_value(self, position1, position2, value):
        self._values[position1][position2] = value

    def add_value(self, position1, position2, value):
        self._safe_append(position1, position2, value)


def load_ABAs(filename):
    ABAs = {}
    ABA_error = {}
    with open(filename) as f:
        line = next(f)
        assert line.startswith('# Target:')
        target = line.strip().split(': ')[1]
        line = next(f)
        assert line.startswith('# Neg Control')
        neg_control_target = line.strip().split(': ')[1]
        line = next(f)
        assert line.startswith('# Concentration')
        line = next(f)
        while not line.startswith('#'):
            max_concentration = float(line.strip().split()[0])
            line = next(f)
        assert line.startswith('# Seq')
        for line in f:
            if line.startswith('#'):
                continue
            words = line.strip().split()
            seq = words[0]
            assert seq not in ABAs, "Duplicate sequence found in ABA file: {}".format(seq)
            ABA, ABA_err = map(float, words[3:])
            ABAs[seq] = max(ABA, 0.0)
            ABA_error[seq] = ABA_err
    if not ABAs:
        print("Warning: no ABAs found!")
    return ABAs, ABA_error


class Comparator(object):
    """
    produces single matrices, already merged, for various types of polymorphisms
    also keep in mind you need to do that scatterplot
    also you need to be able to compare different types of data within a single experiment (e.g. insertions vs mismatches)

    generalize the comparison plots or write all three kinds

    """
    def __init__(self):
        self._experiments = {}

    def add_experiment(self, label, target_sequence, ABAs, ABA_errors):
        self._experiments[label] = {'ts': target_sequence, 'ABAs': ABAs, 'ABA_errors': ABA_errors}

    def compare1d(self, experiment1, experiment2, type1, type2, guide_only=False, normalize=False):
        pass

    def compare_2d_mismatches(self, experiment1, experiment2, guide_only=False, normalize=False, side='lower'):
        return self.compare2d(experiment1, experiment2, 'mismatches', 'mismatches', guide_only=guide_only, normalize=normalize, side=side)

    def compare_2d_insertions(self, experiment1, experiment2, guide_only=False, normalize=False, side='lower'):
        return self.compare2d(experiment1, experiment2, 'insertions', 'insertions', guide_only=guide_only, normalize=normalize, side=side)

    def compare_2d_deletions(self, experiment1, experiment2, guide_only=False, normalize=False, side='lower'):
        return self.compare2d(experiment1, experiment2, 'deletions', 'deletions', guide_only=guide_only, normalize=normalize, side=side)

    def compare_2d_complement_stretches(self, experiment1, experiment2, guide_only=False, normalize=False, side='lower'):
        return self.compare2d(experiment1, experiment2, 'complement_stretches', 'complement_stretches', guide_only=guide_only, normalize=normalize, side=side)

    def compare2d(self, experiment1, experiment2, type1, type2, guide_only=False, normalize=False, return_each_matrix=False, side='lower'):
        """
        This method is mostly used internally, but it does permit the user to compare two different sequence types.
        For example, it would allow you to look at a single experiment and compare mismatches to insertions.

        """
        assert type1 in ('mismatches', 'insertions', 'deletions', 'complement_stretches'), 'Invalid experiment type: %s' % type1
        assert type2 in ('mismatches', 'insertions', 'deletions', 'complement_stretches'), 'Invalid experiment type: %s' % type2

        ABAs1 = self._experiments[experiment1]['ABAs']
        ABAs2 = self._experiments[experiment2]['ABAs']
        ABA_errors1 = self._experiments[experiment1]['ABA_errors']
        ABA_errors2 = self._experiments[experiment2]['ABA_errors']

        if normalize:
            normalize_by1 = ABAs1.get(self._experiments[experiment1]['ts'].sequence) or max(ABAs1.values())
            normalize_by2 = ABAs2.get(self._experiments[experiment2]['ts'].sequence) or max(ABAs2.values())
        else:
            normalize_by1 = None
            normalize_by2 = None

        matrix = self._determine_matrix_type(experiment1, experiment2, type1, type2)
        merge_positions = not self._directly_comparable(experiment1, experiment2, type1, type2)
        flip_sequence = self._experiments[experiment1]['ts'].pam_side != self._experiments[experiment2]['ts'].pam_side

        if guide_only:
            display_sequence1 = self._experiments[experiment1]['ts'].guide.sequence
            display_sequence2 = self._experiments[experiment2]['ts'].guide.sequence
        else:
            display_sequence1 = self._experiments[experiment1]['ts'].sequence
            display_sequence2 = self._experiments[experiment2]['ts'].sequence

        # figure out what sequence to use for plotting and if we can show actual bases or just positions
        return_sequence = display_sequence1 if len(display_sequence1) < len(display_sequence2) else display_sequence2
        if merge_positions:
            sequence_labels = [str(i + 1) for i in range(len(return_sequence))]
        else:
            sequence_labels = ['$%s_{%d}$' % (base, i + 1) for i, base in enumerate(return_sequence)]

        # if one sequence is longer than the other (which will happen if the "sequence" is the same
        # but the PAM is a different length)
        sequence_length = min(len(display_sequence1), len(display_sequence2))
        em1 = matrix(display_sequence1[:sequence_length])
        em2 = matrix(display_sequence2[:sequence_length])
        errm1 = matrix(display_sequence1[:sequence_length])
        errm2 = matrix(display_sequence2[:sequence_length])

        load_func = {'mismatches': self._load_2d_mismatches,
                     'insertions': self._load_2d_insertions,
                     'deletions': self._load_2d_deletions,
                     'complement_stretches': self._load_2d_complement_stretches}

        load_func[type1](em1, errm1, ABAs1, ABA_errors1, self._experiments[experiment1]['ts'], guide_only,
                         sequence_length, merge_positions)
        load_func[type2](em2, errm2, ABAs2, ABA_errors2, self._experiments[experiment2]['ts'], guide_only,
                         sequence_length, merge_positions)
        if not return_each_matrix:
            return return_sequence, sequence_labels, merge_positions, em1.to_matrix(normalize_by=normalize_by1, side=side) - em2.to_matrix(flip_sequence=flip_sequence, normalize_by=normalize_by2, side=side)
        else:
            return em1.to_matrix(normalize_by=normalize_by1, side=side), \
                   em2.to_matrix(flip_sequence=flip_sequence, normalize_by=normalize_by2, side=side), \
                   errm1.to_matrix(normalize_by=normalize_by1, side=side), \
                   errm2.to_matrix(flip_sequence=flip_sequence, normalize_by=normalize_by2, side=side)

    def _load_2d_mismatches(self, matrix, error_matrix, ABAs, ABA_errors, target_sequence, guide_only, sequence_length, merge_positions):
        iterable = target_sequence.guide if guide_only else target_sequence
        for i, j, base_i, base_j, seq in iterable.double_mismatches:
            if i >= sequence_length:
                continue
            if guide_only:
                sequence = target_sequence.pam + seq if target_sequence.pam_side == 5 else seq + target_sequence.pam
            else:
                sequence = seq
            affinity = ABAs.get(sequence)
            error = ABA_errors.get(sequence)
            if merge_positions:
                matrix.add_value(i, j, affinity)
                error_matrix.add_value(i, j, error)
            else:
                matrix.set_value(i, j, base_i, base_j, affinity)
                error_matrix.set_value(i, j, base_i, base_j, error)

    def _load_2d_insertions(self, matrix, error_matrix, ABAs, ABA_errors, target_sequence, guide_only, sequence_length, merge_positions):
        iterable = target_sequence.guide if guide_only else target_sequence
        for i, j, base_i, base_j, seq in iterable.double_insertions:
            if i >= sequence_length:
                continue
            if guide_only:
                sequence = target_sequence.pam + seq if target_sequence.pam_side == 5 else seq + target_sequence.pam
            else:
                sequence = seq
            affinity = ABAs.get(sequence)
            error = ABA_errors.get(sequence)
            if merge_positions:
                matrix.add_value(i, j, affinity)
                error_matrix.add_value(i, j, error)
            else:
                matrix.set_value(i, j, base_i, base_j, affinity)
                error_matrix.set_value(i, j, base_i, base_j, error)

    def _load_2d_deletions(self, matrix, error_matrix, ABAs, ABA_errors, target_sequence, guide_only, sequence_length, merge_positions):
        iterable = target_sequence.guide if guide_only else target_sequence
        for i, j, seq in iterable.double_deletions:
            if i >= sequence_length:
                continue
            if guide_only:
                sequence = target_sequence.pam + seq if target_sequence.pam_side == 5 else seq + target_sequence.pam
            else:
                sequence = seq
            affinity = ABAs.get(sequence)
            error = ABA_errors.get(sequence)
            if merge_positions:
                matrix.add_value(i, j, affinity)
                error_matrix.add_value(i, j, error)
            else:
                matrix.set_value(i, j, affinity)
                error_matrix.set_value(i, j, error)

    def _load_2d_complement_stretches(self, matrix, error_matrix, ABAs, ABA_errors, target_sequence, guide_only, sequence_length, merge_positions):
        iterable = target_sequence.guide if guide_only else target_sequence
        for start, stop, seq in iterable.complement_stretches:
            if start >= sequence_length:
                continue
            if guide_only:
                sequence = target_sequence.pam + seq if target_sequence.pam_side == 5 else seq + target_sequence.pam
            else:
                sequence = seq
            affinity = ABAs.get(sequence)
            error = ABA_errors.get(sequence)
            if merge_positions:
                matrix.add_value(stop, start, affinity)
                error_matrix.add_value(stop, start, error)
            else:
                matrix.set_value(stop, start, affinity)
                error_matrix.set_value(stop, start, error)

    def _directly_comparable(self, experiment1, experiment2, type1, type2):
        if type1 != type2:
            return False
        directions_same = self._experiments[experiment1]['ts'].pam_side == self._experiments[experiment2]['ts'].pam_side
        sequence_same = self._experiments[experiment1]['ts'].sequence == self._experiments[experiment2]['ts'].sequence
        if not directions_same or not sequence_same:
            return False
        return True

    def _determine_matrix_type(self, experiment1, experiment2, type1, type2):
        if not self._directly_comparable(experiment1, experiment2, type1, type2):
            return SinglePositionMatrix
        if type1 in ('complement_stretches', 'deletions'):
            return SinglePositionMatrix
        if type1 == 'mismatches':
            return MismatchMatrix
        if type1 == 'insertions':
            return InsertionMatrix
        raise ValueError("Could not determine matrix type!")


def load_config_value(item_name, override_value):
    # let the user override this method with a manually-specified value
    if override_value:
        return override_value
    try:
        with open("champ.yml") as f:
            config = yaml.load(f)
            return config[item_name]
    except Exception as e:
        print(e)
        raise ValueError("We could not determine the {item_name} from champ.yml. "
                         "Make sure you have a config file and that the value is set.".format(item_name=item_name))


def determine_flow_cell_id(override_value):
    # let the user override this method with a manually-specified value
    if override_value:
        return override_value

    flow_cell_regex = re.compile(r'.*?(SA\d{5}).*?')
    filenames = glob.glob('*')
    candidates = set()
    # look through all the filenames in this directory and look for anything that looks like a sequencing run ID,
    # which we use as flow cell IDs
    for filename in filenames:
        match = flow_cell_regex.search(filename)
        if match:
            candidates.add(match.group(1))
    # there should be a unique value unless someone had an unfortunate choice of filenames
    if len(candidates) == 1:
        return list(candidates)[0]
    raise ValueError(
        "We're unable to automatically determine the flow cell ID based on filenames, you'll need to set it manually.")


def determine_data_channel(all_channels, alignment_channel):
    alignment_channel = set([alignment_channel])
    all_channels = set(map(str, all_channels))
    data_channel = all_channels - alignment_channel
    if len(data_channel) == 1:
        # There are two channels, so we return the one that isn't the alignment channel
        return list(data_channel)[0]
    if not data_channel:
        # There is only one channel, so alignment markers and data are in the same channel
        return list(alignment_channel)[0]
    raise ValueError("Could not determine data channel. You'll need to set this manually.")
