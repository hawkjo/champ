from champ.kd import filter_reads_with_unusual_intensities, assemble_read_intensities_for_fitting
import numpy as np


def test_filter_reads_with_unusual_intensities():
    a = [1, 2, 4, 8, 16, 32, 64]
    b = [1, 2, 4, 8, 16, 32, 64]
    c = [1, 2, 4, 8, 16, 32, 64]
    d = [1, 2, 4, 98, 16, 32, 64]
    e = [1, 2, 4, 8, 16, 32, 64]
    f = [1, 2, 4, 8, 16, 32, 64]
    g = [1, 2, 4, 8, 16, 32, 64]
    intensities = [a, b, c, d, e, f, g]
    good_intensities = filter_reads_with_unusual_intensities(intensities)
    assert len(good_intensities) == 6
    for i in good_intensities:
        assert i[3] == 8


def test_assemble_read_intensities_for_fitting():
    read_name_intensities = {'a': [1,      2,       3,  4,  5],
                             'b': [2,      np.nan,  3,  4, 77],
                             'c': [np.nan, np.nan, 34, 35, 36],
                             'd': [3, 3, 3, 3, 3],
                             'e': [1, 2, 3, 4, 6]}
    intensities = assemble_read_intensities_for_fitting(['a', 'b', 'c'], read_name_intensities)
    assert intensities == [[1, 2], [2], [3, 3, 34], [4, 4, 35], [5, 77, 36]]
