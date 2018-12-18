from champ.constants import MINIMUM_REQUIRED_COUNTS
from collections import defaultdict
import numpy as np
import h5py
from scipy.optimize import curve_fit
from sklearn.neighbors import KernelDensity
from scipy.optimize import minimize
import progressbar
import lomp


BOOTSTRAP_ROUNDS = 30
MAX_BOOTSTRAP_SAMPLE_SIZE = 2000
TUKEY_CONSTANT = 1.5


def fit_kd(x, Kd):
    """ Requres normalized intensities. """
    return 1.0 / (1 + (float(Kd) / x))


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


def get_minimum_intensity(neg_control_intensities):
    background_intensities = [i[0] for i in neg_control_intensities if not np.isnan(i[0])]
    return get_mode(background_intensities)


def get_maximum_intensity(fit_func, concentrations, matched_intensities):
    all_concentrations, all_intensities = assemble_flat_concentrations_and_intensities(concentrations, matched_intensities)
    popt, _ = curve_fit(fit_func, all_concentrations, all_intensities)
    kd, imax = popt
    return kd, imax


def determine_kd(fitting_concentrations, fitting_intensities):
    try:
        popt, pcov = curve_fit(fit_kd, fitting_concentrations, fitting_intensities)
        kd = popt[0]
    except (FloatingPointError, RuntimeError, Exception) as e:
        return None
    else:
        return kd


def bootstrap_kd(all_concentrations, normalized_intensities, return_all_bootstrapped_kds=False, rounds=None, cluster_count=None):
    rounds = BOOTSTRAP_ROUNDS if rounds is None else rounds
    cluster_count = MAX_BOOTSTRAP_SAMPLE_SIZE if cluster_count is None else cluster_count
    kds = []
    # in case some fits don't work, do some extra rounds until we have the number we want
    for i in range(rounds*10):
        sample_of_intensities = sample_lists_with_replacement(normalized_intensities)[:cluster_count]
        fitting_concentrations, fitting_intensities = assemble_flat_concentrations_and_intensities(all_concentrations, sample_of_intensities)
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


def find_boundary_parameters(concentrations, neg_control_intensities, matched_intensities):
    imin = get_minimum_intensity(neg_control_intensities)

    def fit_kd_and_imax(x, Kd, Imax):
        return (Imax - imin) / (1.0 + (float(Kd) / x)) + imin

    matched_kd, imax = get_maximum_intensity(fit_kd_and_imax, concentrations, matched_intensities)
    all_neg_intensities = normalize_intensities(neg_control_intensities, imin, imax)
    neg_conc, neg_int = assemble_flat_concentrations_and_intensities(concentrations, all_neg_intensities)
    (neg_kd,), pcov = curve_fit(fit_kd, neg_conc, neg_int)
    return matched_kd, neg_kd, imin, imax


def _thread_fit_kd(sequence_intensities, all_concentrations, imin, imax):
    # group_intensities is a tuple of a unique label (typically a sequence of interest or location in the genome)
    # and intensities is a list of lists, with each member being the value of an intensity gradient
    sequence, intensities = sequence_intensities
    normalized_intensities = [normalize_intensities(ints, imin, imax) for ints in intensities]
    normalized_intensities = filter_reads_with_insufficient_observations(normalized_intensities, len(all_concentrations) - 3)
    normalized_intensities = filter_reads_with_unusual_intensities(normalized_intensities)
    if len(normalized_intensities) < MINIMUM_REQUIRED_COUNTS:
        return None
    fitting_concentrations, fitting_intensities = assemble_flat_concentrations_and_intensities(all_concentrations, normalized_intensities)
    kd = determine_kd(fitting_concentrations, fitting_intensities)
    if kd is None:
        return None
    kd_uncertainty = bootstrap_kd(all_concentrations, normalized_intensities)
    if kd_uncertainty is None:
        return None
    return sequence, kd, kd_uncertainty, len(intensities)


def fit_all_kds(sequence_read_name_intensities, all_concentrations, imin, imax, process_count=8):
    # sequence_read_name_intensities: List[Dict[str, List[List[float]]]
    # sequence_read_name_intensities should be a list of dictionaries that map read names to intensities
    # each dictionary should all be related to some group of reads that have the same sequence or overlap the same
    # region of the genome
    for result in lomp.parallel_map(sequence_read_name_intensities.items(),
                                    _thread_fit_kd,
                                    args=(all_concentrations, imin, imax),
                                    process_count=process_count):
        if result is not None:
            yield result


def calculate_all_synthetic_kds(h5_filename, concentrations, interesting_reads_filename, matched_sequence, neg_control_sequence):
    string_dt = h5py.special_dtype(vlen=str)
    kd_dt = np.dtype([('sequence', string_dt),
                      ('kd', np.float),
                      ('kd_uncertainty', np.float),
                      ('count', np.int32)])
    print("about to open")
    with h5py.File(h5_filename, 'a') as h5:

        read_name_intensities = load_read_name_intensities(h5_filename)
        interesting_read_names = load_sequence_read_names(interesting_reads_filename)
        sequence_read_name_intensities = defaultdict(list)
        for sequence, read_names in interesting_read_names.items():
            for read_name in read_names:
                if read_name not in read_name_intensities:
                    continue
                sequence_read_name_intensities[sequence].append(read_name_intensities[read_name])

        matched_intensities = filter_reads_with_unusual_intensities(sequence_read_name_intensities[matched_sequence])
        neg_control_intensities = filter_reads_with_unusual_intensities(sequence_read_name_intensities[neg_control_sequence])
        matched_kd, neg_kd, imin, imax = find_boundary_parameters(concentrations, neg_control_intensities, matched_intensities)

        dataset = h5.create_dataset('synthetic-kds', (1,), dtype=kd_dt, maxshape=(None,))
        index = 0

        with progressbar.ProgressBar(max_value=len(sequence_read_name_intensities)) as pbar:
            for sequence, kd, kd_uncertainty, count in pbar(fit_all_kds(sequence_read_name_intensities, concentrations, imin, imax, process_count=8)):
                if count >= MINIMUM_REQUIRED_COUNTS:
                    dataset.resize((index+1,))
                    dataset[index] = (sequence, kd, kd_uncertainty, count)
                    index += 1


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
