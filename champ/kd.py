from scipy.optimize import curve_fit
import numpy as np
import lomp

BOOTSTRAP_ROUNDS = 20
MAX_BOOTSTRAP_SAMPLE_SIZE = 2000
MINIMUM_READ_COUNT = 5
TUKEY_CONSTANT = 1.5


def hyperbola(concentrations, yint, delta_y, kd):
    return ((delta_y - yint) / (1.0 + (kd / concentrations))) + yint


def fixed_delta_y_hyperbola(delta_y):
    """ We need to pass a hyperbolic function to scipy.optimize.curve_fit, but the value of delta_y has to be hard
    coded and not fit by the algorithm. Here, we build a function that has delta_y baked in and then return the
    function."""
    def func(concentrations, yint, kd):
        return ((delta_y - yint) / (1 + (kd / concentrations))) + yint
    return func


def fit_hyperbola(concentrations, signals, delta_y=None):
    """
    :param concentrations: X-axis values representing concentrations in arbitrary units
    :param signals: Y-axis values representing some kind of signal. Don't normalize this.
    Neither of these can be batched - these are flat 1D lists.

    :return:
        yint: the Y-intercept of the fit, often the background signal
        yint_stddev: standard deviation of the error of yint
        delta_y: total height of the fit
        delta_y_stddev: standard deviation of the error of delta_y
        kd: the dissociation constant
        kd_stddev: the standard deviation of the error of kd

    """
    if delta_y is None:
        (yint, fit_delta_y, kd), _ = curve_fit(hyperbola,
                                               concentrations,
                                               signals,
                                               bounds=((0, 0.0, 10**-100),
                                                   (np.inf, np.inf, np.inf)))
    else:
        func = fixed_delta_y_hyperbola(delta_y)
        (yint, kd), _ = curve_fit(func,
                                  concentrations,
                                  signals,
                                  bounds=((0, 10 ** -100), (np.inf, np.inf)))
        fit_delta_y = delta_y
    return yint, fit_delta_y, kd


def fit_all_kds(group_intensities, all_concentrations, process_count=8, delta_y=None):
    # sequence_read_name_intensities: List[Dict[str, List[List[float]]]
    # sequence_read_name_intensities should be a list of dictionaries that map read names to intensities
    # each dictionary should all be related to some group of reads that have the same sequence or overlap the same
    # region of the genome
    minimum_required_observations = max(len(all_concentrations) - 3, 5)
    for result in lomp.parallel_map(group_intensities.items(),
                                    _thread_fit_kd,
                                    args=(all_concentrations, minimum_required_observations, delta_y),
                                    process_count=process_count):
        if result is not None:
            yield result


def fit_one_group_kd(intensities, all_concentrations, delta_y=None):
    minimum_required_observations = max(len(all_concentrations) - 3, 5)
    try:
        result = _thread_fit_kd((None, intensities),
                                all_concentrations,
                                minimum_required_observations,
                                delta_y)
    except Exception as e:
        return None
    else:
        if result is None:
            return None
        _, kd, kd_uncertainty, yint, fit_delta_y, count = result
        return kd, kd_uncertainty, yint, fit_delta_y, count


def _thread_fit_kd(group_intensities, all_concentrations, minimum_required_observations, delta_y):
    # group_intensities is a tuple of a unique label (typically a sequence of interest or location in the genome)
    # and intensities is a list of lists, with each member being the value of an intensity gradient
    group_unique_label, intensities = group_intensities
    intensities = filter_reads_with_unusual_intensities(intensities)
    fitting_concentrations = []
    fitting_intensities = []
    for intensity_gradient in intensities:
        for intensity, concentration in zip(intensity_gradient, all_concentrations):
            if np.isnan(intensity):
                continue
            fitting_intensities.append(intensity)
            fitting_concentrations.append(concentration)
    if len(set(fitting_concentrations)) < minimum_required_observations:
        return None

    kd, yint, fit_delta_y = fit_kd(fitting_concentrations, fitting_intensities, delta_y=delta_y)
    kd_uncertainty = bootstrap_kd_uncertainty(all_concentrations, intensities, delta_y=delta_y)
    if kd is None or kd_uncertainty is None:
        return None
    return group_unique_label, kd, kd_uncertainty, yint, fit_delta_y, len(intensities)


def fit_kd(all_concentrations, all_intensities, delta_y=None):
    """ all_intensities is a list of dicts, with read_name: intensity"""
    try:
        yint, fit_delta_y, kd = fit_hyperbola(all_concentrations, all_intensities, delta_y=delta_y)
    except (FloatingPointError, RuntimeError, Exception) as e:
        return None, None, None
    else:
        return kd, yint, fit_delta_y


def sample_lists_with_replacement(lists):
    # there is no random sampling algorithm with replacement in the standard library, and numpy's
    # random.choice requires 1D arrays
    indexes = np.random.randint(len(lists), size=min(MAX_BOOTSTRAP_SAMPLE_SIZE, len(lists)))
    return [lists[index] for index in indexes]


def bootstrap_kd_uncertainty(all_concentrations, all_intensities, delta_y=None):
    kds = []
    for i in range(BOOTSTRAP_ROUNDS):
        sample_of_intensities = sample_lists_with_replacement(all_intensities)
        intensities = []
        concentrations = []
        for n, concentration in enumerate(list(all_concentrations)):
            for intensity_gradient in sample_of_intensities:
                if n < len(intensity_gradient):
                    intensity = intensity_gradient[n]
                    if np.isnan(intensity):
                        continue
                    intensities.append(intensity)
                    concentrations.append(concentration)
        try:
            _, _, kd = fit_hyperbola(concentrations, intensities, delta_y=delta_y)
        except (FloatingPointError, RuntimeError, Exception) as e:
            continue
        else:
            kds.append(kd)
    if not kds:
        return None
    return np.std(kds)


def filter_reads_with_unusual_intensities(intensities):
    """
    Filters out intensity gradients where any of the measurements were absurdly high or low. Each intensity gradient
    is from a single cluster of DNA.

    :param intensities: a list of numpy arrays, potentially with np.nan values

    """
    bad_indexes = set()
    assert len(set([len(intensity) for intensity in intensities])) == 1, "All reads should have the same number of observations. Missing observations should be represented by np.nan"
    for index in range(len(intensities[0])):
        index_intensities = [intensity_gradient[index] for intensity_gradient in intensities if np.isnan(intensity_gradient[index])]
        if len(index_intensities) < MINIMUM_READ_COUNT:
            continue
        q1 = np.percentile(index_intensities, 25)
        q3 = np.percentile(index_intensities, 75)
        iqr = q3 - q1
        min_range, max_range = (q1 - TUKEY_CONSTANT * iqr, q3 + TUKEY_CONSTANT * iqr)
        for n, intensity_gradient in enumerate(intensities):
            if intensity_gradient[index] is not np.nan and (intensity_gradient[index] < min_range or intensity_gradient[index] > max_range):
                bad_indexes.add(n)
    return [ints for n, ints in enumerate(intensities) if n not in bad_indexes]
