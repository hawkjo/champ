import glob
import os
from collections import defaultdict
import sys

import flabpal
import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import progressbar
import yaml
from scipy import stats

from champ import misc, intensity, initialize, seqtools, interactive, git_commit
from champ.constants import MINIMUM_REQUIRED_COUNTS
from champ.gbtools import genome_main
from champ.kd import fit_all_kds, saturated_at_concentration, fit_one_group_kd
from champ.seqtools import build_interesting_sequences

try:
    matplotlib.style.use('flab')
except:
    pass


process_count = int(sys.argv[1]) if len(sys.argv) > 1 else 8
flow_cell_id = interactive.determine_flow_cell_id()
target_name = interactive.load_config_value('perfect_target_name')
neg_control_target_name = interactive.load_config_value('neg_control_target_name')
all_channels = list(map(str, initialize.determine_channel_names('.')))
alignment_channel = interactive.load_config_value('alignment_channel')
data_channel = interactive.determine_data_channel(all_channels, alignment_channel)
nonneg_lda_weights_fpath = interactive.load_config_value('lda_weights', default='/shared/yeast_beast_LDA_weights.txt')
pam_size = int(interactive.load_config_value('pam_size', default=4))
extended_pam_size = int(interactive.load_config_value('extended_pam_size', default=6))
pam_side = int(interactive.load_config_value('pam_side', default=5))

read_name_dir = os.path.join('/shared', flow_cell_id, 'read_names')
bamfile_path = os.path.join('/shared', flow_cell_id, 'all_fastqs', 'genomic.bam')
read_name_kd_filename = os.path.join('results', 'cluster-data.h5')
read_names_by_seq_fpath = os.path.join(read_name_dir, 'read_names_by_seq.txt')

with open('/shared/targets.yml') as f:
    targets = yaml.load(f)

target = targets[target_name]
neg_control_target = targets[neg_control_target_name]

print('Flow Cell ID: {}'.format(flow_cell_id))
print('Target {}: {}'.format(target_name, target))
print('Neg control target {}: {}'.format(neg_control_target_name, neg_control_target))
print('Channels: {}'.format(all_channels))
print('Protein channel: {}'.format(data_channel))
print('Output file: {}'.format(read_name_kd_filename))
print('Git commit used for this analysis: {}'.format(git_commit))

interesting_seqs = set()

stretch = set()
for i in range(1, len(target) + 1):
    stretch.update(seqtools.get_stretch_of_complement_seqs(target, i))
insertions = set()
for i in range(1, 3):
    insertions.update(seqtools.get_contiguous_insertion_seqs(target, i))
for i in range(1, 3):
    insertions.update(seqtools.get_insertion_seqs(target, i))
deletions = set()
for i in range(1, 3):
    deletions.update(seqtools.get_deletion_seqs(target, i))
mismatches = set()
for i in range(1, 3):
    mismatches.update(seqtools.get_mismatch_seqs(target, i))
pam_end = '5p' if pam_side == 5 else '3p'
extended_pam = seqtools.get_randomized_pam_seqs(target, pam_size, extended_pam_size, end=pam_end)
other_targets = set()
for s in targets.values():
    other_targets.add(s)

interesting_seqs.update(other_targets)
interesting_seqs.update(stretch)
interesting_seqs.update(insertions)
interesting_seqs.update(deletions)
interesting_seqs.update(mismatches)
interesting_seqs.update(extended_pam)

try:
    with open("additional-sequences-of-interest.txt") as f:
        for line in f:
            line = line.strip()
            if line:
                interesting_seqs.add(line)
except IOError:
    pass

print("Interesting sequences: %d" % len(interesting_seqs))

interesting_read_names_filename = os.path.join(read_name_dir, 'interesting_{target_name}_reads_by_seq.txt'.format(
    target_name=target_name.lower()))
if os.path.exists(interesting_read_names_filename):
    # No need to recalculate, we can just load this from disk
    interesting_read_names = {}
    with open(interesting_read_names_filename) as f:
        for line in f:
            line = line.split("\t")
            sequence = line[0]
            read_names = line[1:]
            interesting_read_names[sequence] = read_names
else:
    interesting_read_names = build_interesting_sequences(read_names_by_seq_fpath, interesting_seqs)
    with open(interesting_read_names_filename, 'w') as f:
        for sequence, read_names in interesting_read_names.items():
            f.write("%s\t%s\n" % (sequence, "\t".join(read_names)))

print("Found read names for %d sequences of interest." % len(interesting_read_names))

all_read_name_fpath = os.path.join(read_name_dir, 'all_read_names.txt')
target_read_name_fpath = os.path.join(read_name_dir, 'target_{}_read_names.txt'.format(target_name.lower()))
perfect_target_read_name_fpath = os.path.join(read_name_dir, 'perfect_target_{}_read_names.txt'.format(target_name.lower()))
neg_control_target_read_name_fpath = os.path.join(read_name_dir, 'perfect_target_{}_read_names.txt'.format(neg_control_target_name.lower()))
phiX_read_name_fpath = os.path.join(read_name_dir, 'phix_read_names.txt')

all_read_names = set(line.strip() for line in open(all_read_name_fpath))
print("All read names: %d" % len(all_read_names))
target_read_names = set(line.strip() for line in open(target_read_name_fpath))
print("Target read names: %d" % len(target_read_names))
perfect_target_read_names = set(line.strip() for line in open(perfect_target_read_name_fpath))
print("Perfect target read names: %d" % len(perfect_target_read_names))
neg_control_target_read_names = set(line.strip() for line in open(neg_control_target_read_name_fpath))
print("Negative control read names: %d" % len(neg_control_target_read_names))
phiX_read_names = set(line.strip() for line in open(phiX_read_name_fpath))
print("Phix read names: %d" % len(phiX_read_names))

h5_fpaths = glob.glob('*.h5')
h5_fpaths.sort(key=misc.parse_concentration)
for fpath in h5_fpaths:
    print misc.parse_concentration(fpath), fpath

results_dirs = [
    os.path.join('results',
                 os.path.splitext(os.path.basename(h5_fpath))[0])
    for h5_fpath in h5_fpaths
]
for d in results_dirs:
    print(d)

int_scores = intensity.IntensityScores(h5_fpaths)
int_scores.get_LDA_scores(results_dirs, nonneg_lda_weights_fpath)
int_scores.normalize_scores()
int_scores.plot_aligned_images('br', 'o*')
int_scores.plot_normalization_constants()
int_scores.print_reads_per_channel()

# The number of observations we require to consider a cluster valid enough for fitting.
# Clusters at the edge of a field of view might not be visible in every concentration due to
# random imperfections in the motion of the stage, and some fields of view might simply not
# align under each concentration.
minimum_observations_per_cluster = len(h5_fpaths) - 3

int_scores.build_good_read_names(minimum_observations_per_cluster)
good_read_names = int_scores.good_read_names
good_perfect_read_names = perfect_target_read_names & good_read_names
print('Good Reads:', len(good_read_names))
print('Good Perfect Reads:', len(good_perfect_read_names))

int_scores.build_score_given_read_name_given_channel()

aligned_read_names = set()
for h5_fpath in h5_fpaths:
    for d in int_scores.scores[h5_fpath][data_channel].values():
        for read_name in d.keys():
            aligned_read_names.add(read_name)
print('Aligned reads in data channel: {}'.format(len(aligned_read_names)))

read_name_intensities = defaultdict(list)
for h5_fpath in h5_fpaths:
    score_given_read_name = int_scores.score_given_read_name_in_channel[h5_fpath][data_channel]
    for read_name in aligned_read_names:
        intensity_val = score_given_read_name.get(read_name, np.nan)
        read_name_intensities[read_name].append(intensity_val)

intensity_matrix = []
# now that aligned_read_names are guaranteed to be unique, we need to give them a guaranteed
# order to that they can stay aligned with the intensity_matrix
aligned_read_names = list(aligned_read_names)
for read_name in aligned_read_names:
    intensity_matrix.append(read_name_intensities[read_name])

intensity_matrix = np.array(intensity_matrix)
string_dt = h5py.special_dtype(vlen=str)

with h5py.File(read_name_kd_filename, 'w') as h5:
    h5.create_dataset('read_names', data=aligned_read_names, dtype=string_dt)

with h5py.File(read_name_kd_filename, 'a') as h5:
    h5.create_dataset('intensities', data=intensity_matrix)

# Find the KD of the perfect target so we can determine at what concentrations the clusters
# should be saturated
sequence_read_name_intensities = defaultdict(list)
for sequence, read_names in interesting_read_names.items():
    for read_name in read_names:
        if read_name not in read_name_intensities:
            continue
        sequence_read_name_intensities[sequence].append(read_name_intensities[read_name])

all_concentrations = [misc.parse_concentration(h5_fpath) for h5_fpath in h5_fpaths]
perfect_kd, perfect_kd_uncertainty, perfect_yint, perfect_delta_y, perfect_counts = fit_one_group_kd(sequence_read_name_intensities[target], all_concentrations, delta_y=None)
print("Perfect target KD is %.1f +/- %.3f nM" % (perfect_kd, perfect_kd_uncertainty))

neg_kd, neg_kd_uncertainty, neg_yint, neg_delta_y, neg_counts = fit_one_group_kd(sequence_read_name_intensities[targets['D']], all_concentrations, delta_y=None)
print("Neg target KD is %.1f +/- %.3f nM" % (neg_kd, neg_kd_uncertainty))

# Determine the median intensity of a saturated cluster
saturated = saturated_at_concentration(perfect_kd)
print("Should be saturated at %.1f nM" % saturated)
saturated_indexes = [index for index, concentration in enumerate(all_concentrations) if concentration > saturated]
if not saturated_indexes:
    # this should never happen, but we'll try to take something reasonable
    print("Warning: perfect target sequence probably did not saturate its target!")
    saturated_indexes = [len(all_concentrations) - 1]

saturated_intensities = []
for intensity_gradient in sequence_read_name_intensities[target]:
    for index in saturated_indexes:
        try:
            value = intensity_gradient[index]
            if not np.isnan(value):
                saturated_intensities.append(value)
        except IndexError:
            continue
median_saturated_intensity = np.median(saturated_intensities)
print("Median saturated intensity: %d (N=%d)" % (median_saturated_intensity, len(saturated_intensities)))

fig, ax = plt.subplots(figsize=(8,8))
ax.plot(all_concentrations, sequence_read_name_intensities[target][0], color=flabpal.blue, alpha=0.005, label='Perfect target intensities')
for intensity_gradient in sequence_read_name_intensities[target][1:2000]:
    ax.plot(all_concentrations, intensity_gradient, color=flabpal.blue, alpha=0.005)
ax.axhline(median_saturated_intensity, linestyle='--', color=flabpal.red, label='Median saturated intensity')
ax.set_title("Perfect target intensities")
ax.set_ylabel("Intensity (AU)")
ax.set_xlabel("Concentration (nM)")
legend = ax.legend()
for handle in legend.legendHandles:
    handle.set_alpha(1.0)

string_dt = h5py.special_dtype(vlen=str)
kd_dt = np.dtype([('sequence', string_dt),
                  ('kd', np.float),
                  ('kd_uncertainty', np.float),
                  ('y_intercept', np.float),
                  ('delta_y', np.float),
                  ('count', np.int32)])
with h5py.File(read_name_kd_filename, 'a') as h5:
    dataset = h5.create_dataset('synthetic-kds', (1,), dtype=kd_dt, maxshape=(None,))
    index = 0
    with progressbar.ProgressBar(max_value=len(sequence_read_name_intensities)) as pbar:
        for sequence, kd, kd_uncertainty, yint, delta_y, count in pbar(fit_all_kds(sequence_read_name_intensities, all_concentrations, process_count=process_count, delta_y=median_saturated_intensity)):
            if count >= MINIMUM_REQUIRED_COUNTS:
                dataset.resize((index+1,))
                dataset[index] = (sequence, kd, kd_uncertainty, yint, delta_y, count)
                index += 1

print("Determined KDs for %d of %d sequences" % (index, len(sequence_read_name_intensities)))

fig, (kd_ax, count_ax) = plt.subplots(1, 2, figsize=(16, 6))

with h5py.File(read_name_kd_filename, 'r') as h5:
    histogram_kds = [i[1] for i in h5['synthetic-kds'][:]]
    kd_ax.hist(histogram_kds, bins=100)
    kd_ax.set_ylabel("Count")
    kd_ax.set_xlabel("$K_D$ (nM)")
    kd_ax.set_title("Histogram of Synthetic Target KDs")

with h5py.File(read_name_kd_filename, 'r') as h5:
    histogram_counts = [i[5] for i in h5['synthetic-kds'][:]]
    count_median = np.median(histogram_counts)
    count_iqr = stats.iqr(histogram_counts)
    count_ax.hist(histogram_counts, bins=20, range=(0, int(count_median + count_iqr * 5)))
    count_ax.set_ylabel("Count")
    count_ax.set_xlabel("Clusters per Sequence")
    count_ax.set_title("Representation of Synthetic Sequences")
fig.tight_layout()

genome_main(bamfile_path, read_name_kd_filename, all_concentrations, median_saturated_intensity, process_count=process_count)

with h5py.File(read_name_kd_filename, 'a') as h5:
    group = h5.create_group("git-commit")
    group.create_group(git_commit)
