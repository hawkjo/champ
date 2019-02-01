from champ.constants import MINIMUM_REQUIRED_COUNTS
from collections import defaultdict
import numpy as np
import h5py
from scipy.optimize import curve_fit, minimize
from sklearn.neighbors import KernelDensity
import progressbar
import lomp
import random


BOOTSTRAP_ROUNDS = 100
MAX_BOOTSTRAP_SAMPLE_SIZE = 2000
TUKEY_CONSTANT = 1.5


def calculate_all_synthetic_kds(h5_filename, concentrations, interesting_read_names, matched_sequence,
                                neg_control_sequence, process_count):
    read_name_intensities = load_read_name_intensities(h5_filename)
    sequence_read_name_intensities = defaultdict(list)
    for sequence, read_names in interesting_read_names.items():
        for read_name in read_names:
            if read_name not in read_name_intensities:
                continue
            sequence_read_name_intensities[sequence].append(read_name_intensities[read_name])
    perfect_kd, perfect_kd_uncertainty, perfect_yint, perfect_delta_y, perfect_counts = fit_one_group_kd(
        sequence_read_name_intensities[matched_sequence], concentrations, delta_y=None)
    print("Perfect target KD is %.1f +/- %.3f nM" % (perfect_kd, perfect_kd_uncertainty))

    neg_kd, neg_kd_uncertainty, neg_yint, neg_delta_y, neg_counts = fit_one_group_kd(
        sequence_read_name_intensities[neg_control_sequence], concentrations, delta_y=perfect_delta_y)
    print("Neg target KD is %.1f +/- %.3f nM" % (neg_kd, neg_kd_uncertainty))

    # Determine the median intensity of a saturated cluster
    saturated = saturated_at_concentration(perfect_kd)
    print("Should be saturated at %.1f nM" % saturated)
    saturated_indexes = [index for index, concentration in enumerate(concentrations) if
                         concentration > saturated]
    if not saturated_indexes:
        # this should never happen, but we'll try to take something reasonable
        print("Warning: perfect target sequence probably did not saturate its target!")
        saturated_indexes = [len(concentrations) - 1]

    saturated_intensities = []
    for intensity_gradient in sequence_read_name_intensities[matched_sequence]:
        for index in saturated_indexes:
            try:
                value = intensity_gradient[index]
                if not np.isnan(value):
                    saturated_intensities.append(value)
            except IndexError:
                continue
    median_saturated_intensity = int(np.median(saturated_intensities))
    print("Median saturated intensity: %d (N=%d)" % (median_saturated_intensity, len(saturated_intensities)))

    string_dt = h5py.special_dtype(vlen=str)
    kd_dt = np.dtype([('sequence', string_dt),
                      ('kd', np.float),
                      ('kd_uncertainty', np.float),
                      ('count', np.int32)])

    with h5py.File(h5_filename, 'a') as h5:
        dataset = h5.create_dataset('synthetic-kds', (1,), dtype=kd_dt, maxshape=(None,))
        index = 0
        with progressbar.ProgressBar(max_value=len(sequence_read_name_intensities)) as pbar:
            for sequence, kd, kd_uncertainty, yint, delta_y, count in pbar(
                    fit_all_kds(sequence_read_name_intensities, concentrations, process_count=process_count,
                                delta_y=median_saturated_intensity)):
                if count >= MINIMUM_REQUIRED_COUNTS:
                    dataset.resize((index + 1,))
                    dataset[index] = (sequence, kd, kd_uncertainty, count)
                    index += 1
    print("Fit %d sequences" % index)


def delta_aba(kd, perfect_kd):
    ratio = kd / perfect_kd
    return np.log(ratio)


def convert_kd_to_normalized_delta_aba(kd, perfect_kd, neg_kd):
    neg_aba = delta_aba(neg_kd, perfect_kd)
    daba = delta_aba(kd, perfect_kd)
    return daba / neg_aba


def load_sequence_read_names(filename):
    sequence_read_names = {}
    with open(filename) as f:
        for line in f:
            line = line.strip().split('\t')
            sequence = line[0]
            read_names = line[1:]
            sequence_read_names[sequence] = read_names
    return sequence_read_names


def load_read_name_intensities(h5_filename):
    with h5py.File(h5_filename, 'r') as h5:
        intensities = h5['intensities'][:]
        read_names = h5['read_names'][:]
        read_name_intensities = {}
        for n, (read_name, intensity_curve) in enumerate(zip(read_names, intensities)):
            read_name_intensities[read_name] = intensity_curve
    # we should have lots of reads
    return read_name_intensities


def create_sequence_read_name_intensities(read_name_intensities, read_names_by_sequence):
    sequence_read_name_intensities = defaultdict(list)
    for sequence, read_names in read_names_by_sequence.items():
        for read_name in read_names:
            if read_name not in read_name_intensities:
                continue
            sequence_read_name_intensities[sequence].append(read_name_intensities[read_name])
    return sequence_read_name_intensities


def get_mode(vals):
    vals = np.array(vals)
    h = 1.06 * np.std(vals) * len(vals)**(-1.0/5.0)
    kdf = KernelDensity(bandwidth=h)
    kdf.fit(vals.reshape(len(vals), 1))

    def neg_kdf(x):
        return -kdf.score(np.array((x,)))

    res = minimize(neg_kdf, x0=np.median(vals), method='Nelder-Mead')
    assert res.success, res
    return float(res.x)


def assemble_flat_concentrations_and_intensities(all_concentrations, cluster_intensities):
    """ all_concentrations is just a sorted list of the unique concentrations used in an experiment.
    it should have like 7 to 11 total values normally

    cluster_intensities is a list of lists, with each member being an intensity curve for one cluster"""
    concentrations = []
    intensities = []
    for n, concentration in enumerate(list(all_concentrations)):
        for intensity_gradient in cluster_intensities:
            if n < len(intensity_gradient):
                intensity = intensity_gradient[n]
                if np.isnan(intensity):
                    continue
                intensities.append(intensity)
                concentrations.append(concentration)
    return concentrations, intensities


def normalize_intensities(intensities, imin, imax):
    d = imax - imin
    return [(np.array(intensity) - imin) / d for intensity in intensities]


def determine_kd(fitting_concentrations, fitting_intensities):
    try:
        popt, pcov = curve_fit(fit_kd, fitting_concentrations, fitting_intensities)
        kd = popt[0]
    except (FloatingPointError, RuntimeError, Exception) as e:
        return None
    else:
        return kd


def bootstrap_kd(concentrations, normalized_intensities, return_all_bootstrapped_kds=False, rounds=None, cluster_count=None):
    if len(normalized_intensities) == 0:
        return None
    rounds = BOOTSTRAP_ROUNDS if rounds is None else rounds
    cluster_count = MAX_BOOTSTRAP_SAMPLE_SIZE if cluster_count is None else cluster_count
    kds = []
    # in case some fits don't work, do some extra rounds until we have the number we want
    for i in range(rounds*10):
        sample_of_intensities = [random.choice(normalized_intensities) for _ in normalized_intensities][:cluster_count]
        # sample_of_intensities = sample_lists_with_replacement(normalized_intensities)[:cluster_count]
        fitting_concentrations, fitting_intensities = assemble_flat_concentrations_and_intensities(concentrations, sample_of_intensities)
        try:
            popt, pcov = curve_fit(fit_kd, fitting_concentrations, fitting_intensities)
            kd = popt[0]
        except (FloatingPointError, RuntimeError, Exception) as e:
            continue
        else:
            kds.append(kd)
            if len(kds) == rounds:
                if return_all_bootstrapped_kds:
                    return kds
                return np.std(kds)
    return None


def convert_2d_list_to_column_list(values):
    """
    Takes a list of lists, all of equal length, and splits the data by each column.

    """
    return np.array([column for column in np.array(values).T])


def get_quality_normalized_intensities(intensities, concentrations, imin, imax):
    normalized_intensities = [normalize_intensities(ints, imin, imax) for ints in intensities]
    normalized_intensities = filter_reads_with_insufficient_observations(normalized_intensities, len(concentrations) - 3)
    return filter_reads_with_unusual_intensities(normalized_intensities)


def fit_kd(all_concentrations, all_intensities, delta_y=None):
    """ all_intensities is a list of dicts, with read_name: intensity"""
    try:
        yint, fit_delta_y, kd = fit_hyperbola(all_concentrations, all_intensities, delta_y=delta_y)
    except (FloatingPointError, RuntimeError, Exception) as e:
        return None, None, None
    else:
        return kd, yint, fit_delta_y


def fit_one_kd(normalized_intensities, concentrations):
    fitting_concentrations, fitting_intensities = assemble_flat_concentrations_and_intensities(concentrations,
                                                                                               normalized_intensities)
    return determine_kd(fitting_concentrations, fitting_intensities)


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


def _thread_fit_kd(group_intensities, all_concentrations, minimum_required_observations, delta_y, bootstrap=True):
    # group_intensities is a tuple of a unique label (typically a sequence of interest or location in the genome)
    # and intensities is a list of lists, with each member being the value of an intensity gradient
    group_unique_label, intensities = group_intensities
    intensities = filter_reads_with_insufficient_observations(intensities, minimum_required_observations)
    if len(intensities) < MINIMUM_REQUIRED_COUNTS:
        return None
    fitting_concentrations = []
    fitting_intensities = []
    for intensity_gradient in intensities:
        for n, (intensity, concentration) in enumerate(zip(intensity_gradient, all_concentrations)):
            if np.isnan(intensity):
                continue
            fitting_intensities.append(intensity)
            fitting_concentrations.append(concentration)

    kd, yint, fit_delta_y = fit_kd(fitting_concentrations, fitting_intensities, delta_y=delta_y)
    if bootstrap:
        kd_uncertainty = bootstrap_kd_uncertainty(all_concentrations, intensities, delta_y=delta_y)
    else:
        kd_uncertainty = 0.0
    if kd is None or kd_uncertainty is None:
        return None
    return group_unique_label, kd, kd_uncertainty, yint, fit_delta_y, len(intensities)


def fit_one_group_kd(intensities, all_concentrations, delta_y=None, bootstrap=True):
    minimum_required_observations = max(len(all_concentrations) - 3, 5)
    try:
        result = _thread_fit_kd((None, intensities),
                                all_concentrations,
                                minimum_required_observations,
                                delta_y, bootstrap=bootstrap)
    except Exception as e:
        return None
    else:
        if result is None:
            return None
        _, kd, kd_uncertainty, yint, fit_delta_y, count = result
        return kd, kd_uncertainty, yint, fit_delta_y, count


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


def hyperbola(concentrations, yint, delta_y, kd):
    return ((delta_y - yint) / (1.0 + (kd / concentrations))) + yint


def fixed_delta_y_hyperbola(delta_y):
    """ We need to pass a hyperbolic function to scipy.optimize.curve_fit, but the value of delta_y has to be hard
    coded and not fit by the algorithm. Here, we build a function that has delta_y baked in and then return the
    function."""
    def func(concentrations, yint, kd):
        return ((delta_y - yint) / (1 + (kd / concentrations))) + yint
    return func


def saturated_at_concentration(kd):
    """ Determines what concentration should have (mostly) saturated clusters given a KD. """
    saturated_fraction = 0.95
    return float(kd * saturated_fraction)/(1.0 - saturated_fraction)


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
                                               bounds=((0.0, 0.0, 10 ** -280),
                                                   (np.inf, np.inf, np.inf)))
    else:
        func = fixed_delta_y_hyperbola(delta_y)
        (yint, kd), _ = curve_fit(func,
                                  concentrations,
                                  signals,
                                  bounds=((0.0, 10 ** -10), (np.inf, np.inf)))
        fit_delta_y = delta_y
    return yint, fit_delta_y, kd


def filter_reads_with_unusual_intensities(intensities):
    """
    Filters out intensity gradients where any of the measurements were absurdly high or low. Each intensity gradient
    is from a single cluster of DNA.

    :param intensities: a list of numpy arrays, potentially with np.nan values

    """
    if len(intensities) == 0:
        return []
    bad_clusters = set()
    assert len(set([len(intensity) for intensity in intensities])) == 1, "All reads should have the same number of " \
                                                                         "observations. Missing observations should be " \
                                                                         "represented by np.nan. %s" % intensities
    for index in range(len(intensities[0])):
        index_intensities = [intensity_gradient[index] for intensity_gradient in intensities if not np.isnan(intensity_gradient[index])]
        if not index_intensities:
            # all values were np.nan, so we can't use this concentration at all
            continue
        q1 = np.percentile(index_intensities, 25)
        q3 = np.percentile(index_intensities, 75)
        iqr = q3 - q1
        min_range, max_range = (q1 - TUKEY_CONSTANT * iqr, q3 + TUKEY_CONSTANT * iqr)
        for n, intensity_gradient in enumerate(intensities):
            if intensity_gradient[index] is not np.nan and (intensity_gradient[index] < min_range or intensity_gradient[index] > max_range):
                bad_clusters.add(n)
    return [ints for n, ints in enumerate(intensities) if n not in bad_clusters]


def filter_reads_with_insufficient_observations(intensities, minimum_required_observations):
    good_intensities = []
    for gradient in intensities:
        if np.sum(~np.isnan(gradient)) >= minimum_required_observations:
            good_intensities.append(gradient)
    return good_intensities


def sample_lists_with_replacement(lists):
    # there is no random sampling algorithm with replacement in the standard library, and numpy's
    # random.choice requires 1D arrays
    indexes = np.random.randint(len(lists), size=min(MAX_BOOTSTRAP_SAMPLE_SIZE, len(lists)))
    return [lists[index] for index in indexes]


def copy_over_everything_but_kds(h5_filename, new_h5_filename):
    with h5py.File(h5_filename, 'r') as h5:
        intensities = h5['intensities'][:]
        read_names = h5['read_names'][:]

    with h5py.File(new_h5_filename, 'w') as h5:
        h5['intensities'] = intensities
        h5['read_names'] = read_names
