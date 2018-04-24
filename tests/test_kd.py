from champ.kd import filter_reads_with_unusual_intensities, bootstrap_kd_uncertainty, fit_all_kds, hyperbola
import numpy as np


def test_fit_all_kds_with_delta_y():
    concentrations = np.array([0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
    a = hyperbola(concentrations, 100, 120000, 4.0)
    b = hyperbola(concentrations, 200, 120000, 5.0)
    c = hyperbola(concentrations, 100, 120000, 4.2)
    d = hyperbola(concentrations, 200, 120000, 4.5)
    e = hyperbola(concentrations, 300, 120000, 5.5)
    f = hyperbola(concentrations, 200, 120000, 4.7)
    g = hyperbola(concentrations, 300, 120000, 5.7)
    h = hyperbola(concentrations, 200, 120000, 5.3)
    sequence_1_read_intensities = [a, b, c, d]
    sequence_2_read_intensities = [e, f, g, h]
    read_name_intensities = {'sequence1': sequence_1_read_intensities, 'sequence2': sequence_2_read_intensities}
    results = list(fit_all_kds(read_name_intensities, concentrations, delta_y=180000, process_count=1))
    assert len(results) == 2
    for read_name, kd, kd_uncertainty, yint, delta_y in results:
        assert kd > 70.0
        assert 0.0 < kd_uncertainty < 1.0


def test_fit_all_kds():
    concentrations = np.array([0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
    a = hyperbola(concentrations, 100, 120000, 4.0)
    b = hyperbola(concentrations, 200, 120000, 5.0)
    c = hyperbola(concentrations, 100, 120000, 4.2)
    d = hyperbola(concentrations, 200, 120000, 4.5)
    e = hyperbola(concentrations, 300, 120000, 5.5)
    f = hyperbola(concentrations, 200, 120000, 4.7)
    g = hyperbola(concentrations, 300, 120000, 5.7)
    h = hyperbola(concentrations, 200, 120000, 5.3)
    i = hyperbola(concentrations, 800, 120000, 5.9)
    j = hyperbola(concentrations, 20, 120000, 4.1)
    read_name_intensities = {'a': [a, b, c, d, e, f, g, h],
                             'b': [j, i, h, g, f, e, d, c]}
    results = list(fit_all_kds(read_name_intensities, concentrations, process_count=1))
    assert len(results) == 2
    for read_name, kd, kd_uncertainty, yint, delta_y in results:
        assert 10.0 > kd > 3.0
        assert 0.0 < kd_uncertainty < 1.0


def test_bootstrap_kd_uncertainty():
    # This really just tests that the function runs and returns a float
    # Since we're randomly sampling, it's possible that we'll get the same list for every sample, and thus have a
    # standard deviation of 0.0
    concentrations = np.array([0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
    a = hyperbola(concentrations, 0, 1, 4.0)
    b = hyperbola(concentrations, 0, 1, 5.0)
    c = hyperbola(concentrations, 0, 1, 3.0)
    d = hyperbola(concentrations, 0, 1, 4.5)
    std = bootstrap_kd_uncertainty(concentrations, map(list, [a, b, c, d]))
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
