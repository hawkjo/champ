# import matplotlib
# matplotlib.use('agg')
# import sys
# import os
# import matplotlib.pyplot as plt
# import fastqimagealigner
# import local_config
# from matplotlib.lines import Line2D
#
#
# def process_image(image_relpath):
#     #----------------------------------------------------------------------
#     # Process file and dir names
#     #----------------------------------------------------------------------
#     project_name = 'SA15151'
#     image_fpath = os.path.realpath(image_relpath)
#     image_fname = os.path.basename(image_fpath)
#     assert image_fname.endswith('.tif'), image_fname
#     image_bname = os.path.splitext(image_fname)[0]
#     image_dir = os.path.split(image_fpath)[0]  # full path of dir
#     image_dirname = os.path.basename(image_dir)  # only name of dir
#
#     results_dir = os.path.join(local_config.results_dir, '20150725', image_dirname)
#     fig_dir = os.path.join(local_config.fig_dir, '20150725', image_dirname)
#     for d in [results_dir, fig_dir]:
#         if not os.path.isdir(d):
#             os.makedirs(d)
#
#     #----------------------------------------------------------------------
#     # Alignment settings
#     #----------------------------------------------------------------------
#     sexcat_fpath = os.path.join(image_dir, image_bname + '.cat')
#     possible_tile_keys = ['lane1tile2110']  # All images in this tile
#     rotation_est = 180.0
#     fq_w_est = 930.0
#
#     read_name_dir = os.path.join(local_config.fourier_data_dir, project_name, 'read_names')
#     our_control_read_name_fpath = os.path.join(read_name_dir, 'control_read_names.txt')
#     our_active_read_name_fpath = os.path.join(read_name_dir, 'dnazyme_read_names.txt')
#     anywhere_perfect_control_read_name_fpath = os.path.join(read_name_dir,
#                                                             'anywhere_perfect_control_read_names.txt')
#     perfect_active_read_name_fpath = os.path.join(read_name_dir, 'perfect_active_read_names.txt')
#     all_read_name_fpath = os.path.join(read_name_dir, 'all_read_names.txt')
#
#     #----------------------------------------------------------------------
#     # Align
#     #----------------------------------------------------------------------
#     fic = fastqimagealigner.FastqImageAligner(project_name)
#     fic.load_reads(tile_data=local_config.fastq_tiles_given_read_name_fpath(our_control_read_name_fpath))
#     fic.set_image_data(image_fpath, 60, median_normalize=True)
#     fic.set_sexcat_from_file(sexcat_fpath)
#     fic.align(possible_tile_keys, rotation_est, fq_w_est, snr_thresh=1.2, min_hits=15,
#               hit_type=['exclusive', 'good_mutual'])
#
#     #----------------------------------------------------------------------
#     # Copy alignment to different sets of points and get aligned rcs
#     #----------------------------------------------------------------------
#     all_fic = fastqimagealigner.FastqImageAligner(project_name)
#     all_tile_data = local_config.fastq_tiles_given_read_name_fpath(all_read_name_fpath)
#     all_fic.all_reads_fic_from_aligned_fic(fic, tile_data=all_tile_data)
#     all_aligned_rcs = all_fic.hitting_tiles[0].aligned_rcs
#
#     our_control_fic = fastqimagealigner.FastqImageAligner(project_name)
#     our_control_tile_data = local_config.fastq_tiles_given_read_name_fpath(our_control_read_name_fpath)
#     our_control_fic.all_reads_fic_from_aligned_fic(fic, tile_data=our_control_tile_data)
#     our_control_aligned_rcs = our_control_fic.hitting_tiles[0].aligned_rcs
#
#     our_active_fic = fastqimagealigner.FastqImageAligner(project_name)
#     our_active_tile_data = local_config.fastq_tiles_given_read_name_fpath(our_active_read_name_fpath)
#     our_active_fic.all_reads_fic_from_aligned_fic(fic, tile_data=our_active_tile_data)
#     our_active_aligned_rcs = our_active_fic.hitting_tiles[0].aligned_rcs
#
#     perfect_control_fic = fastqimagealigner.FastqImageAligner(project_name)
#     perfect_control_tile_data = local_config.fastq_tiles_given_read_name_fpath(
#         anywhere_perfect_control_read_name_fpath
#     )
#     perfect_control_fic.all_reads_fic_from_aligned_fic(fic, tile_data=perfect_control_tile_data)
#     perfect_control_aligned_rcs = perfect_control_fic.hitting_tiles[0].aligned_rcs
#
#     perfect_active_fic = fastqimagealigner.FastqImageAligner(project_name)
#     perfect_active_tile_data = local_config.fastq_tiles_given_read_name_fpath(perfect_active_read_name_fpath)
#     perfect_active_fic.all_reads_fic_from_aligned_fic(fic, tile_data=perfect_active_tile_data)
#     perfect_active_aligned_rcs = perfect_active_fic.hitting_tiles[0].aligned_rcs
#
#     #----------------------------------------------------------------------
#     # Make images
#     #----------------------------------------------------------------------
#     all_hits_fig_fpath = os.path.join(fig_dir, '%s_all_hits.png' % image_bname)
#     our_reads_fig_fpath = os.path.join(fig_dir, '%s_our_reads.png' % image_bname)
#     perfect_reads_fig_fpath = os.path.join(fig_dir, '%s_perfect_reads.png' % image_bname)
#
#     ax = fic.plot_all_hits()
#     ax.get_figure().savefig(all_hits_fig_fpath)
#
#     im = fic.image_data.im
#
#     def overlay_fig(control_rcs, active_rcs, out_fpath, title):
#         fig, ax = plt.subplots(figsize=(12, 12))
#         alpha = 0.3
#         markersize = 10
#         ax.matshow(im, cmap=plt.get_cmap('Blues'))
#         ax.plot(control_rcs[:, 1], control_rcs[:, 0],
#                 'o', color='red', markersize=markersize, alpha=alpha)
#         ax.plot(active_rcs[:, 1], active_rcs[:, 0],
#                 'o', color='darkgoldenrod', markersize=markersize, alpha=alpha)
#         ax.set_xlim((0, im.shape[0]))
#         ax.set_ylim((im.shape[1], 0))
#         ax.set_title(title)
#
#         control_line = Line2D([], [], color='red', alpha=alpha, marker='o',
#                               markersize=markersize, label='Control')
#         active_line = Line2D([], [], color='darkgoldenrod', alpha=alpha, marker='o',
#                              markersize=markersize, label='DNAzyme')
#         handles = [control_line, active_line]
#         legend = ax.legend(handles=handles)
#         legend.get_frame().set_color('white')
#
#         fig.savefig(out_fpath)
#
#     overlay_fig(our_control_aligned_rcs,
#                 our_active_aligned_rcs,
#                 our_reads_fig_fpath,
#                 '%s Our Reads' % image_bname)
#     overlay_fig(perfect_control_aligned_rcs,
#                 perfect_active_aligned_rcs,
#                 perfect_reads_fig_fpath,
#                 '%s Perfect Reads' % image_bname)
#
#     #----------------------------------------------------------------------
#     # Write point alignments
#     #----------------------------------------------------------------------
#     stats_fpath = os.path.join(results_dir, '%s_stats.txt' % image_bname)
#     fic.write_alignment_stats(stats_fpath)
#
#     all_read_rcs_fpath = os.path.join(results_dir, '%s_all_read_rcs.txt' % image_bname)
#     all_fic.write_read_names_rcs(all_read_rcs_fpath)
#
# if __name__ == '__main__':
#     usage_fmt = '%s <image>' % sys.argv[0]
#     if len(sys.argv) != len(usage_fmt.split()):
#         sys.exit('Usage: ' + usage_fmt)
#     image_relpath = sys.argv[1]
#     process_image(image_relpath)
