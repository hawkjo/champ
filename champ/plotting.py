import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from matplotlib import gridspec
import matplotlib as mpl
import flabpal


def plot_hit_hists(fia, ax=None):
    """ Plots histograms of the different quality cluster alignment categories. """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    non_mut_dists = fia.hit_dists(fia.non_mutual_hits)
    bins = np.linspace(0, max(non_mut_dists), 50)

    if non_mut_dists:
        ax.hist(non_mut_dists, bins, label='Non-mutual hits', normed=True, histtype='step')
    if fia.bad_mutual_hits:
        ax.hist(fia.hit_dists(fia.bad_mutual_hits), bins,
                label='Bad mutual hits', normed=True, histtype='step')
    if fia.good_mutual_hits:
        ax.hist(fia.hit_dists(fia.good_mutual_hits), bins, label='Good mutual hits',
                normed=True, histtype='step')
    if fia.exclusive_hits:
        ax.hist(fia.hit_dists(fia.exclusive_hits), bins, label='Exclusive hits',
                normed=True, histtype='step')
    ax.legend()
    ax.set_title('%s Nearest Neighbor Distance Distributions' % fia.image_data.fname)
    return ax


def plot_hits(fia, hits, color, ax, kwargs={}):
    """ Draws circles where hits are located. """
    for i, j in hits:
        ax.plot([fia.clusters.point_rcs[i, 1], fia.aligned_rcs_in_frame[j, 1]],
                [fia.clusters.point_rcs[i, 0], fia.aligned_rcs_in_frame[j, 0]],
                color=color, **kwargs)
    return ax


def plot_ellipses(fia, ax, alpha=1.0, color=(1, 0, 0)):
    ells = [Ellipse(xy=(pt.c, pt.r), width=3, height=3, angle=0.0)
            for pt in fia.clusters.points]
    for e in ells:
        ax.add_artist(e)
        e.set_alpha(alpha)
        e.set_facecolor(color)


def plot_all_hits(fia, im_kwargs={}, line_kwargs={}, fqpt_kwargs={}, sext_kwargs={},
                 title_kwargs={}, legend_kwargs={}):
    """ 
    Creates a plot of a field of view with the raw microscope image in the background and 
    hit locations drawn over them. Provides a very obvious measure of whether the alignment worked or not. 
    
    """
    fig, ax = plt.subplots(figsize=(15, 15))

    kwargs = {'cmap': plt.get_cmap('Blues')}
    kwargs.update(im_kwargs)
    ax.matshow(fia.image_data.image, **kwargs)

    kwargs = {'color': 'k', 'alpha': 0.3, 'linestyle': '', 'marker': 'o', 'markersize': 3}
    kwargs.update(fqpt_kwargs)
    ax.plot(fia.aligned_rcs_in_frame[:, 1], fia.aligned_rcs_in_frame[:, 0], **kwargs)

    kwargs = {'alpha': 0.6, 'color': 'darkgoldenrod'}
    kwargs.update(sext_kwargs)
    plot_ellipses(fia, ax, **kwargs)

    plot_hits(fia, fia.non_mutual_hits, 'grey', ax, line_kwargs)
    plot_hits(fia, fia.bad_mutual_hits, 'b', ax, line_kwargs)
    plot_hits(fia, fia.good_mutual_hits, 'magenta', ax, line_kwargs)
    plot_hits(fia, fia.exclusive_hits, 'r', ax, line_kwargs)
    ax.set_title('All Hits: %s vs. %s\nRot: %s deg, Fq width: %s um, Scale: %s px/fqu, Corr: %s, SNR: %s'
            % (fia.image_data.fname,
               ','.join(tile.key for tile in fia.hitting_tiles),
               ','.join('%.2f' % tile.rotation_degrees for tile in fia.hitting_tiles),
               ','.join('%.2f' % tile.width for tile in fia.hitting_tiles),
               ','.join('%.5f' % tile.scale for tile in fia.hitting_tiles),
               ','.join('%.1f' % tile.best_max_corr if hasattr(tile, 'best_max_corr') else '0.0' for tile in fia.hitting_tiles),
               ','.join('%.2f' % tile.snr if hasattr(tile, 'snr') else '-' for tile in fia.hitting_tiles),
               ), **title_kwargs)
    ax.set_xlim([0, fia.image_data.image.shape[1]])
    ax.set_ylim([fia.image_data.image.shape[0], 0])

    grey_line = Line2D([], [], color='grey',
            label='Non-mutual hits: %d' % (len(fia.non_mutual_hits)))
    blue_line = Line2D([], [], color='blue',
            label='Bad mutual hits: %d' % (len(fia.bad_mutual_hits)))
    magenta_line = Line2D([], [], color='magenta',
            label='Good mutual hits: %d' % (len(fia.good_mutual_hits)))
    red_line = Line2D([], [], color='red',
            label='Exclusive hits: %d' % (len(fia.exclusive_hits)))
    sexcat_line = Line2D([], [], color='darkgoldenrod', alpha=0.6, marker='o', markersize=10,
            label='Sextractor Ellipses: %d' % (len(fia.clusters.point_rcs)))
    fastq_line = Line2D([], [], color='k', alpha=0.3, marker='o', markersize=10,
            label='Fastq Points: %d' % (len(fia.aligned_rcs_in_frame)))
    handles = [grey_line, blue_line, magenta_line, red_line, sexcat_line, fastq_line]
    legend = ax.legend(handles=handles, **legend_kwargs)
    legend.get_frame().set_color('white')
    return ax


def get_cluster_counts(ia, seq):
    """
    Each concentration could have different numbers of clusters, since they won't all successfully align each time.
    We take the lowest number of clusters in a single concentration and use that as the count.
    This assumes that all images that align in that concentration align in all the others,
    and that the images were taken in the same place on the chip during each acquisition. All lists in the intensity array are the same length
    but contain None when that cluster was not found in a particular concentration.

    """
    cluster_counts = []
    for lol in ia.intensity_lol_given_seq[seq]:
        cluster_counts.append(len([i for i in lol if i is not None]))
    return min(cluster_counts) if cluster_counts else 0


def configure_position_penalty_axes(target, fig, penalty_axes, count_axes, xticklabels, fontsize, tick_fontsize, title,
                                    yaxis_type, base_color, target_name, legend=True):
    if yaxis_type == 'kd':
        yaxis_label = '$K_{d} (nM)$'
    elif yaxis_type == 'ddG':
        yaxis_label = '$\Delta \Delta G\ (K_{B}T)$'
    elif yaxis_type == 'ABA':
        yaxis_label = '$ABA\ (K_{B}T)$'
    else:
        yaxis_label = '?????'

    penalty_axes.xaxis.grid(False)
    penalty_axes.set_xlim((-0.5, len(target)-0.5))
    penalty_axes.set_xticks(range(len(target)))
    penalty_axes.set_xticklabels(xticklabels, fontsize=tick_fontsize)
    count_axes.set_yscale('log')
    count_axes.set_xlim((-0.5, len(target)-0.5))
    count_axes.set_xticks(range(len(target)))
    count_axes.set_xticklabels(xticklabels, fontsize=tick_fontsize)
    ylim = penalty_axes.get_ylim()
    for i, c in enumerate(target):
        # color the background with the correct base
        penalty_axes.fill_between([i-0.5, i+0.5], [ylim[0]]*2, [ylim[1]]*2, color=base_color[c], alpha=0.14)
    penalty_axes.set_ylim(ylim)
    penalty_axes.set_title(title, fontsize=fontsize)
    count_axes.set_title("Unique Clusters Per Mismatch Sequence", fontsize=fontsize)
    count_axes.set_ylabel("Count", fontsize=fontsize)
    count_axes.xaxis.set_ticks_position('none')
    penalty_axes.set_xlabel('Target {target_name} Reference Sequence (Background Color)'.format(target_name=target_name), fontsize=fontsize)
    penalty_axes.set_ylabel(yaxis_label, fontsize=fontsize)
    if legend:
        penalty_axes.legend(loc='best')
    penalty_axes.xaxis.set_ticks_position('none')
    fig.tight_layout()


def build_base_colorcode_axis(ax, sequence, base_colors, vertical=False):
    bases = 'ACGT'
    colors = [(1, 1, 1)] + [base_colors[base] for base in bases]
    color_index = {base: i for i, base in enumerate(' ' + bases)}
    cmap = mpl.colors.ListedColormap(colors)
    if not vertical:
        base_data = np.array([[color_index[base] for base in sequence]])
    else:
        data = []
        for base in sequence:
            data.append([color_index[base]])
        base_data = np.array(data)

    mat = ax.matshow(base_data, cmap=cmap, vmin=0, vmax=4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    return ax


def sum_nan_arrays(a, b):
    # used to combine an upper and lower triangle matrix. If the matrices have their diagonal values they will sum (not good)!
    ma = np.isnan(a)
    mb = np.isnan(b)
    return np.where(ma & mb, np.nan, np.where(ma, 0, a) + np.where(mb, 0, b))


def plot_2d_mismatches(sequence, human_readable_indexes, protein_name, lower_ABA_matrix, base_colors=None, upper_ABA_matrix=None, fontsize=18):
    if base_colors is None:
        base_colors = {'A': flabpal.blue, 'C': flabpal.yellow, 'G': flabpal.green, 'T': flabpal.red}
    dimension = 3
    width_ratios = [.5, 1, len(sequence) * dimension, 3]
    height_ratios = [len(sequence) * dimension, 1, .5]
    fig = plt.figure(figsize=(sum(width_ratios) / 3, sum(height_ratios) / 3))

    gs = gridspec.GridSpec(3, 4,
                           width_ratios=width_ratios,
                           height_ratios=height_ratios,
                           wspace=0.01, hspace=0.01
                           )
    data_index = 2
    left_sequence_index = 0
    bottom_sequence_index = 10
    left_color_index = 1
    bottom_color_index = 6
    cbar_index = 3

    sequence_labels = ["$%s_{%d}$" % (base, index) for base, index in human_readable_indexes]

    # Add the sequence labels to the left of the figure
    left_sequence_ax = fig.add_subplot(gs[left_sequence_index])
    left_sequence_ax.set_yticklabels(sequence_labels[::-1], fontsize=8*dimension)
    left_sequence_ax.set_yticks([dimension * x + dimension / 2.0 for x in range(len(sequence))])
    left_sequence_ax.set_ylim([0, len(sequence) * dimension])
    left_sequence_ax.spines['top'].set_visible(False)
    left_sequence_ax.spines['right'].set_visible(False)
    left_sequence_ax.spines['bottom'].set_visible(False)
    left_sequence_ax.spines['left'].set_visible(False)
    left_sequence_ax.tick_params(top="off")
    left_sequence_ax.tick_params(bottom="off")
    left_sequence_ax.tick_params(right="off")
    left_sequence_ax.tick_params(left="off")
    left_sequence_ax.set_xticklabels([])

    # Add the sequence labels to the bottom of the figure
    bottom_sequence_ax = fig.add_subplot(gs[bottom_sequence_index])
    bottom_sequence_ax.set_xticklabels(sequence_labels, fontsize=8*dimension)
    bottom_sequence_ax.set_xticks([dimension * x + 1.5 for x in range(len(sequence))])
    bottom_sequence_ax.set_xlim([0, len(sequence) * dimension - 1])
    bottom_sequence_ax.spines['top'].set_visible(False)
    bottom_sequence_ax.spines['right'].set_visible(False)
    bottom_sequence_ax.spines['bottom'].set_visible(False)
    bottom_sequence_ax.spines['left'].set_visible(False)
    bottom_sequence_ax.tick_params(top="off")
    bottom_sequence_ax.tick_params(bottom="off")
    bottom_sequence_ax.tick_params(right="off")
    bottom_sequence_ax.tick_params(left="off")
    bottom_sequence_ax.set_yticklabels([])

    # Add the color bars to the left and bottom of the figure to indicate which base the mismatch has been converted to
    mm_bases = ''.join(['ACGT'.replace(base, '') for base in sequence])
    left_color_codes_ax = fig.add_subplot(gs[left_color_index])
    build_base_colorcode_axis(left_color_codes_ax, mm_bases, base_colors, vertical=True)
    bottom_color_codes_ax = fig.add_subplot(gs[bottom_color_index])
    build_base_colorcode_axis(bottom_color_codes_ax, mm_bases, base_colors)

    # Add data to the main part of the figure
    data_ax = fig.add_subplot(gs[data_index])
    data_ax.set_axis_bgcolor(0.87 * np.array([1, 1, 1]))
    if upper_ABA_matrix is None:
        ms = data_ax.matshow(lower_ABA_matrix, cmap='viridis')
    else:
        # we "add" the arrays, retaining NaNs, to create a comparison matrix
        # both matrices should have their include_diagonal_values=False or else those will sum,
        # or if one is on it will be misleading
        ms = data_ax.matshow(sum_nan_arrays(upper_ABA_matrix, lower_ABA_matrix), cmap='viridis')
    data_ax.set_yticks([])
    data_ax.set_xticks([])
    fig.suptitle('%s Double Mismatch Binding Affinities' % protein_name, fontsize=fontsize*3)

    # Add a color bar to the right side to quantify the colors in the main figure
    cbar_ax = fig.add_subplot(gs[cbar_index])
    cbar_ax.tick_params(labelsize=18)
    cbar = plt.colorbar(ms, cax=cbar_ax)
    cbar.set_label('$ABA (k_{B}T)$', fontsize=fontsize*2)


def plot_2d_deletions(sequence, sequence_labels, lower_ABA_matrix, upper_ABA_matrix=None, fontsize=18):
    width_ratios = [.5, len(sequence), 1]
    height_ratios = [len(sequence), .5]
    fig = plt.figure(figsize=(sum(width_ratios), sum(height_ratios)))

    gs = gridspec.GridSpec(2, 3,
                           width_ratios=width_ratios,
                           height_ratios=height_ratios,
                           wspace=0.01, hspace=0.01)
    data_index = 1
    left_sequence_index = 0
    bottom_sequence_index = 4
    cbar_index = 2

    # Add the sequence labels to the left of the figure
    left_sequence_ax = fig.add_subplot(gs[left_sequence_index])
    left_sequence_ax.set_yticklabels(sequence_labels[::-1], fontsize=fontsize)
    left_sequence_ax.set_yticks([x + 0.5 for x in range(len(sequence))])
    left_sequence_ax.set_ylim([0, len(sequence)])
    left_sequence_ax.spines['top'].set_visible(False)
    left_sequence_ax.spines['right'].set_visible(False)
    left_sequence_ax.spines['bottom'].set_visible(False)
    left_sequence_ax.spines['left'].set_visible(False)
    left_sequence_ax.tick_params(top="off")
    left_sequence_ax.tick_params(bottom="off")
    left_sequence_ax.tick_params(right="off")
    left_sequence_ax.tick_params(left="off")
    left_sequence_ax.set_xticklabels([])

    # Add the sequence labels to the bottom of the figure
    bottom_sequence_ax = fig.add_subplot(gs[bottom_sequence_index])
    bottom_sequence_ax.set_xticklabels(sequence_labels, fontsize=fontsize)
    bottom_sequence_ax.set_xticks([x + 0.5 for x in range(len(sequence))])
    bottom_sequence_ax.set_xlim([0, len(sequence)])
    bottom_sequence_ax.spines['top'].set_visible(False)
    bottom_sequence_ax.spines['right'].set_visible(False)
    bottom_sequence_ax.spines['bottom'].set_visible(False)
    bottom_sequence_ax.spines['left'].set_visible(False)
    bottom_sequence_ax.tick_params(top="off")
    bottom_sequence_ax.tick_params(bottom="off")
    bottom_sequence_ax.tick_params(right="off")
    bottom_sequence_ax.tick_params(left="off")
    bottom_sequence_ax.set_yticklabels([])

    # Add data to the main part of the figure
    data_ax = fig.add_subplot(gs[data_index])
    data_ax.set_axis_bgcolor(0.87 * np.array([1, 1, 1]))
    if upper_ABA_matrix is None:
        ms = data_ax.matshow(lower_ABA_matrix, cmap='viridis')
    else:
        # we "add" the arrays, retaining NaNs, to create a comparison matrix
        # both matrices should have their include_diagonal_values=False or else those will sum,
        # and if one is on it will be misleading
        ms = data_ax.matshow(sum_nan_arrays(upper_ABA_matrix, lower_ABA_matrix), cmap='viridis')
    data_ax.set_yticks([])
    data_ax.set_xticks([])

    # Add a color bar to the right side to quantify the colors in the main figure
    cbar_ax = fig.add_subplot(gs[cbar_index])
    cbar_ax.tick_params(labelsize=18)
    cbar = plt.colorbar(ms, cax=cbar_ax)
    cbar.set_label('$ABA (k_{B}T)$', fontsize=fontsize*2)
