import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse


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
