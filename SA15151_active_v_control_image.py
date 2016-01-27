import matplotlib
matplotlib.use('agg')
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import random
import tifffile
import fastqimagecorrelator
import local_config
from matplotlib.lines import Line2D
import nd2reader


def process_image(nd2_fpath, im_idx):
    #----------------------------------------------------------------------
    # Process file and dir names
    #----------------------------------------------------------------------
    project_name = 'SA15151'
    vmax=1

    nd2_fname = os.path.basename(nd2_fpath)
    assert nd2_fname.endswith('.nd2'), nd2_fname
    nd2_bname = os.path.splitext(nd2_fname)[0]
    nd2_parent_dir = os.path.split(nd2_fpath)[0]  # full path of dir
    date = os.path.split(nd2_parent_dir)[-1]
    assert date.startswith('20') and len(date) == 8 and set(date) <= set('0123456789'), date

    results_dir = os.path.join(local_config.results_dir, date, nd2_fname)
    fig_dir = os.path.join(local_config.fig_dir, date, nd2_fname)
    for d in [results_dir, fig_dir]:
        if not os.path.isdir(d):
            os.makedirs(d)

    #----------------------------------------------------------------------
    # Alignment settings
    #----------------------------------------------------------------------
    nd2 = nd2reader.Nd2(nd2_fpath)
    sexcat_fpath = os.path.join(nd2_parent_dir, nd2_bname, '%d.cat' % im_idx)
    possible_tile_keys = ['lane1tile%d' % x for x in range(2109, 2114)]
    rotation_est = 180.0
    fq_w_est = 930.0

    read_name_dir = os.path.join(local_config.fourier_data_dir, project_name, 'read_names')
    our_control_read_name_fpath = os.path.join(read_name_dir, 'control_read_names.txt')
    our_active_read_name_fpath = os.path.join(read_name_dir, 'dnazyme_read_names.txt')
    anywhere_perfect_control_read_name_fpath = os.path.join(read_name_dir,
                                                            'anywhere_perfect_control_read_names.txt')
    perfect_active_read_name_fpath = os.path.join(read_name_dir, 'perfect_active_read_names.txt')
    all_read_name_fpath = os.path.join(read_name_dir, 'all_read_names.txt')

    #----------------------------------------------------------------------
    # Align
    #----------------------------------------------------------------------
    fic = fastqimagecorrelator.FastqImageAligner(project_name)
    fic.load_reads(tile_data=local_config.fastq_tiles_given_read_name_fpath(our_control_read_name_fpath))
    fic.set_image_data(im=nd2[im_idx].data, objective=60, median_normalize=True)
    fic.set_sexcat_from_file(sexcat_fpath)
    fic.align(possible_tile_keys, rotation_est, fq_w_est, snr_thresh=1.2, min_hits=15,
              hit_type=['exclusive', 'good_mutual'])

    #----------------------------------------------------------------------
    # Copy alignment to different sets of points and get aligned rcs
    #----------------------------------------------------------------------
    all_fic = fastqimagecorrelator.FastqImageAligner(project_name)
    all_tile_data = local_config.fastq_tiles_given_read_name_fpath(all_read_name_fpath)
    all_fic.all_reads_fic_from_aligned_fic(fic, tile_data=all_tile_data)
    all_aligned_rcs = all_fic.hitting_tiles[0].aligned_rcs

    our_control_fic = fastqimagecorrelator.FastqImageAligner(project_name)
    our_control_tile_data = local_config.fastq_tiles_given_read_name_fpath(our_control_read_name_fpath)
    our_control_fic.all_reads_fic_from_aligned_fic(fic, tile_data=our_control_tile_data)
    our_control_aligned_rcs = our_control_fic.hitting_tiles[0].aligned_rcs

    our_active_fic = fastqimagecorrelator.FastqImageAligner(project_name)
    our_active_tile_data = local_config.fastq_tiles_given_read_name_fpath(our_active_read_name_fpath)
    our_active_fic.all_reads_fic_from_aligned_fic(fic, tile_data=our_active_tile_data)
    our_active_aligned_rcs = our_active_fic.hitting_tiles[0].aligned_rcs

    perfect_control_fic = fastqimagecorrelator.FastqImageAligner(project_name)
    perfect_control_tile_data = local_config.fastq_tiles_given_read_name_fpath(
        anywhere_perfect_control_read_name_fpath
    )
    perfect_control_fic.all_reads_fic_from_aligned_fic(fic, tile_data=perfect_control_tile_data)
    perfect_control_aligned_rcs = perfect_control_fic.hitting_tiles[0].aligned_rcs

    perfect_active_fic = fastqimagecorrelator.FastqImageAligner(project_name)
    perfect_active_tile_data = local_config.fastq_tiles_given_read_name_fpath(perfect_active_read_name_fpath)
    perfect_active_fic.all_reads_fic_from_aligned_fic(fic, tile_data=perfect_active_tile_data)
    perfect_active_aligned_rcs = perfect_active_fic.hitting_tiles[0].aligned_rcs

    #----------------------------------------------------------------------
    # Make images
    #----------------------------------------------------------------------
    all_hits_fig_fpath = os.path.join(fig_dir, '%s_all_hits.png' % im_idx)
    our_reads_fig_fpath = os.path.join(fig_dir, '%s_our_reads.png' % im_idx)
    perfect_reads_fig_fpath = os.path.join(fig_dir, '%s_perfect_reads.png' % im_idx)

    ax = fic.plot_all_hits()
    ax.get_figure().savefig(all_hits_fig_fpath)

    im = fic.image_data.im

    def overlay_fig(control_rcs, active_rcs, out_fpath, title):
        fig, ax = plt.subplots(figsize=(12, 12))
        alpha = 0.3
        markersize = 10
        ax.matshow(im, cmap=plt.get_cmap('Blues'), vmax=vmax)
        ax.plot(control_rcs[:, 1], control_rcs[:, 0],
                'o', color='red', markersize=markersize, alpha=alpha)
        ax.plot(active_rcs[:, 1], active_rcs[:, 0],
                'o', color='darkgoldenrod', markersize=markersize, alpha=alpha)
        ax.set_xlim((0, im.shape[0]))
        ax.set_ylim((im.shape[1], 0))
        ax.set_title(title)

        control_line = Line2D([], [], color='red', alpha=alpha, marker='o',
                              markersize=markersize, label='Control')
        active_line = Line2D([], [], color='darkgoldenrod', alpha=alpha, marker='o',
                             markersize=markersize, label='DNAzyme')
        handles = [control_line, active_line]
        legend = ax.legend(handles=handles)
        legend.get_frame().set_color('white')

        fig.savefig(out_fpath)

    overlay_fig(our_control_aligned_rcs,
                our_active_aligned_rcs,
                our_reads_fig_fpath,
                '%s Our Reads' % im_idx)
    overlay_fig(perfect_control_aligned_rcs,
                perfect_active_aligned_rcs,
                perfect_reads_fig_fpath,
                '%s Perfect Reads' % im_idx)
    
    #----------------------------------------------------------------------
    # Write point alignments
    #----------------------------------------------------------------------
    stats_fpath = os.path.join(results_dir, '%s_stats.txt' % im_idx)
    fic.write_alignment_stats(stats_fpath)

    all_read_rcs_fpath = os.path.join(results_dir, '%s_all_read_rcs.txt' % im_idx)
    all_fic.write_read_names_rcs(all_read_rcs_fpath)

if __name__ == '__main__':
    usage_fmt = '%s <nd2_fpath> <image_idx>' % sys.argv[0]
    if len(sys.argv) != len(usage_fmt.split()):
        sys.exit('Usage: ' + usage_fmt)
    nd2_fpath = sys.argv[1]
    im_idx = int(sys.argv[2])
    process_image(nd2_fpath, im_idx)
