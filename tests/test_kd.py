from champ.kd import filter_reads_with_unusual_intensities, assemble_read_intensities_for_fitting, \
    assemble_fitting_inputs, bootstrap_kd_uncertainty
import numpy as np
from biofits import hyperbola


def test_bootstrap_kd_uncertainty():
    # This really just tests that the function runs and returns a float
    # Since we're randomly sampling, it's possible that we'll get the same list for every sample, and thus have a
    # standard deviation of 0.0
    concentrations = np.array([0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
    a = hyperbola(concentrations, 0, 1, 4.0)
    b = hyperbola(concentrations, 0, 1, 5.0)
    c = hyperbola(concentrations, 0, 1, 3.0)
    d = hyperbola(concentrations, 0, 1, 4.5)
    std = bootstrap_kd_uncertainty(concentrations, [a, b, c, d])
    assert std >= 0.0


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


def test_filter_reads_with_unusual_intensities_some_nans():
    a = [1, 2, 4, 8, 16, 32, np.nan]
    b = [1, 2, 4, 8, 16, 32, 64]
    c = [1, np.nan, 4, 8, 16, 32, 64]
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
    read_name_intensities = [[1,      2,       3,  4,  5],
                             [2,      np.nan,  3,  4, 77],
                             [np.nan, np.nan, 34, 35, 36]]
    intensities = assemble_read_intensities_for_fitting(read_name_intensities)
    assert intensities == [[1, 2], [2], [3, 3, 34], [4, 4, 35], [5, 77, 36]]


def test_assemble_fitting_inputs():
    all_concentrations = [1, 2, 4, 8, 16]
    assembled_intensities = [[3, 4, 3, 2], [], [9, 5, 6, 7], [21, 22, 21, 20, 19], []]
    concentrations, concentrations_per_observation, intensities = assemble_fitting_inputs(assembled_intensities,
                                                                                          all_concentrations)
    assert concentrations == [1, 4, 8]
    assert intensities == [[3, 4, 3, 2], [9, 5, 6, 7], [21, 22, 21, 20, 19]]
    assert concentrations_per_observation == [1, 1, 1, 1, 4, 4, 4, 4, 8, 8, 8, 8, 8]
