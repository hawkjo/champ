import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse


def plot_hit_hists(fastq_image_aligner, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    non_mut_dists = fastq_image_aligner.hit_dists(fastq_image_aligner.non_mutual_hits)
    bins = np.linspace(0, max(non_mut_dists), 50)

    if non_mut_dists:
        ax.hist(non_mut_dists, bins, label='Non-mutual hits', normed=True, histtype='step')
    if fastq_image_aligner.bad_mutual_hits:
        ax.hist(fastq_image_aligner.hit_dists(fastq_image_aligner.bad_mutual_hits), bins, label='Bad mutual hits', normed=True, histtype='step')
    if fastq_image_aligner.good_mutual_hits:
        ax.hist(fastq_image_aligner.hit_dists(fastq_image_aligner.good_mutual_hits), bins, label='Good mutual hits', normed=True, histtype='step')
    if fastq_image_aligner.exclusive_hits:
        ax.hist(fastq_image_aligner.hit_dists(fastq_image_aligner.exclusive_hits), bins, label='Exclusive hits', normed=True, histtype='step')
    ax.legend()
    ax.set_title('%s Nearest Neighbor Distance Distributions' % fastq_image_aligner.image_data.bname)
    return ax


def plot_threshold_gmm(fastq_image_aligner, axs=None, force=False):
    if axs is None:
        fig, axs = plt.subplots(1, 2, figsize=(15, 6))
    non_mut_dists = fastq_image_aligner.hit_dists(fastq_image_aligner.non_mutual_hits)
    if not hasattr(fastq_image_aligner, 'gmm') and force:
        fastq_image_aligner.gmm_thresh(non_mut_dists)
    xs = np.linspace(0, max(non_mut_dists), 200)
    posteriors = fastq_image_aligner.gmm.predict_proba(xs)
    pdf = np.exp(fastq_image_aligner.gmm.score_samples(xs)[0])

    axs[0].hist(non_mut_dists, 40, histtype='step', normed=True, label='Data')
    axs[0].plot(xs, pdf, label='PDF')
    ylim = axs[0].get_ylim()
    axs[0].plot([fastq_image_aligner.second_neighbor_thresh, fastq_image_aligner.second_neighbor_thresh], ylim, 'g--', label='Threshold')
    axs[0].set_title('%s GMM PDF of Non-mutual hits' % fastq_image_aligner.image_data.bname)
    axs[0].legend()
    axs[0].set_ylim(ylim)

    axs[1].hist(non_mut_dists, 40, histtype='step', normed=True, label='Data')
    axs[1].plot(xs, posteriors, label='Posterior')
    axs[1].plot([fastq_image_aligner.second_neighbor_thresh, fastq_image_aligner.second_neighbor_thresh], [0, 1], 'g--', label='Threshold')
    axs[1].set_title('%s GMM Posterior Probabilities' % fastq_image_aligner.image_data.bname)
    axs[1].legend()
    return axs


def plot_hits(fastq_image_aligner, hits, color, ax, kwargs):
    for i, j in hits:
        ax.plot([fastq_image_aligner.sexcat.point_rcs[i, 1], fastq_image_aligner.aligned_rcs_in_frame[j, 1]],
                [fastq_image_aligner.sexcat.point_rcs[i, 0], fastq_image_aligner.aligned_rcs_in_frame[j, 0]],
                color=color, **kwargs)
    return ax


def plot_sextraction_ellipses(sextraction, ax=None, alpha=1.0, color=(1, 0, 0)):
        if ax is None:
            fig, ax = plt.subplots()
        ells = [Ellipse(xy=(pt.c, pt.r), width=pt.width, height=pt.height, angle=pt.theta) for pt in self.points]
        for e in ells:
            ax.add_artist(e)
            e.set_alpha(alpha)
            e.set_facecolor(color)


def plot_all_hits(fastq_image_aligner, ax=None, im_kwargs=None, line_kwargs=None, fqpt_kwargs=None, sext_kwargs=None,
                 title_kwargs=None, legend_kwargs=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(15, 15))

    kwargs = {'cmap': plt.get_cmap('Blues')}
    kwargs.update(im_kwargs or {})
    ax.matshow(fastq_image_aligner.image_data.im, **kwargs)

    kwargs = {'color': 'k', 'alpha': 0.3, 'linestyle': '', 'marker': 'o', 'markersize': 3}
    kwargs.update(fqpt_kwargs or {})
    ax.plot(fastq_image_aligner.aligned_rcs_in_frame[:, 1], fastq_image_aligner.aligned_rcs_in_frame[:, 0], **kwargs)

    kwargs = {'alpha': 0.6, 'color': 'darkgoldenrod'}
    kwargs.update(sext_kwargs or {})
    plot_sextraction_ellipses(fastq_image_aligner.sexcat, ax=ax, **kwargs)

    line_kwargs = line_kwargs or {}
    title_kwargs = title_kwargs or {}
    legend_kwargs = legend_kwargs or {}
    fastq_image_aligner.plot_hits(fastq_image_aligner.non_mutual_hits, 'grey', ax, line_kwargs)
    fastq_image_aligner.plot_hits(fastq_image_aligner.bad_mutual_hits, 'b', ax, line_kwargs)
    fastq_image_aligner.plot_hits(fastq_image_aligner.good_mutual_hits, 'magenta', ax, line_kwargs)
    fastq_image_aligner.plot_hits(fastq_image_aligner.exclusive_hits, 'r', ax, line_kwargs)
    ax.set_title('All Hits: %s vs. %s %s\nRot: %s deg, Fq width: %s um, Scale: %s px/fqu, Corr: %s, SNR: %s'
            % (fastq_image_aligner.image_data.bname,
               fastq_image_aligner.chip_id,
               ','.join(tile.key for tile in fastq_image_aligner.hitting_tiles),
               ','.join('%.2f' % tile.rotation_degrees for tile in fastq_image_aligner.hitting_tiles),
               ','.join('%.2f' % tile.w for tile in fastq_image_aligner.hitting_tiles),
               ','.join('%.5f' % tile.scale for tile in fastq_image_aligner.hitting_tiles),
               ','.join('%.1f' % tile.best_max_corr for tile in fastq_image_aligner.hitting_tiles),
               ','.join('%.2f' % tile.snr if hasattr(tile, 'snr') else '-' for tile in fastq_image_aligner.hitting_tiles),
               ), **title_kwargs)
    ax.set_xlim([0, fastq_image_aligner.image_data.im.shape[1]])
    ax.set_ylim([fastq_image_aligner.image_data.im.shape[0], 0])

    grey_line = Line2D([], [], color='grey', label='Non-mutual hits: %d' % (len(fastq_image_aligner.non_mutual_hits)))
    blue_line = Line2D([], [], color='blue', label='Bad mutual hits: %d' % (len(fastq_image_aligner.bad_mutual_hits)))
    magenta_line = Line2D([], [], color='magenta', label='Good mutual hits: %d' % (len(fastq_image_aligner.good_mutual_hits)))
    red_line = Line2D([], [], color='red', label='Exclusive hits: %d' % (len(fastq_image_aligner.exclusive_hits)))
    sexcat_line = Line2D([], [], color='darkgoldenrod', alpha=0.6, marker='o', markersize=10,
                         label='Sextractor Ellipses: %d' % (len(fastq_image_aligner.sexcat.point_rcs)))
    fastq_line = Line2D([], [], color='k', alpha=0.3, marker='o', markersize=10,
                        label='Fastq Points: %d' % (len(fastq_image_aligner.aligned_rcs_in_frame)))
    handles = [grey_line, blue_line, magenta_line, red_line, sexcat_line, fastq_line]
    legend = ax.legend(handles=handles, **legend_kwargs)
    legend.get_frame().set_color('white')
    return ax


def plot_hit_vectors(fastq_image_aligner, hit_types=('exclusive',), ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(15, 15))
    colors = {'exclusive': 'r',
              'good_mutual': 'magenta',
              'bad_mutual': 'b',
              'non_mutual': 'grey'}
    for hit_type in hit_types:
        hits = getattr(fastq_image_aligner, hit_type + '_hits')
        pts = np.array([fastq_image_aligner.sexcat.point_rcs[i] - fastq_image_aligner.aligned_rcs_in_frame[j] for i, j in hits])
        ax.plot(pts[:, 1], pts[:, 0], '.', color=colors[hit_type])
    ax.plot([0], [0], 'k*')
    ax.set_aspect(1)
    ylim = ax.get_ylim()
    ax.set_ylim((ylim[1], ylim[0]))
    ax.set_title('{0} {1} Hit Diffs'.format(fastq_image_aligner.image_data.bname, hit_type.capitalize()))
    ax.set_xlabel('c')
    ax.set_ylabel('r')
    return ax
