from champ.kd import filter_reads_with_unusual_intensities


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
