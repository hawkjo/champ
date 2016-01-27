# import matplotlib
# matplotlib.use('agg')
# import sys
# import numpy as np
# import matplotlib.pyplot as plt
# import os
# import glob
# import random
# import tifffile
# from matplotlib import animation
# from JSAnimation import IPython_display
# import fastqimagecorrelator
# import local_config
# import nd2reader
# import string
#
#
# date = '20150816'
# datadir = os.path.join(local_config.data_dir, date)
# project_name = 'SA15161'
#
# nd2_fpath = '/home/hawkjo/IlyaProjects/miseq_alignment/data/20150816/20150816-SA15161_blockingP7_SP1_attoshort_AlexaSTV_pepQD706.nd2'
# nd2 = nd2reader.Nd2(nd2_fpath)
#
# fq_dir = os.path.join(local_config.data_dir, 'from_fourierseq', project_name, 'all_fastqs')
# all_read_name_fpath = os.path.join(fq_dir, 'all_read_names.txt')
# sexcat_dir = os.path.splitext(nd2._filename)[0]
# possible_tile_keys = ['lane1tile%d' % i for i in range(2104, 2108)]
# rotation_est = 180
# fq_w_est = 935
#
# locations = [[0, (445.0, 465.0, 200.0, 220.0)],
#               [16, (325.0, 345.0, 385.0, 405.0)],
#               [22, (80.0, 100.0, 240.0, 260.0)],
#               [56, (285.0, 305.0, 490.0, 510.0)],
#               [62,
#                  (235.0, 255.0, 25.0, 45.0),
#                  (395.0, 415.0, 55.0, 75.0),
#                  (270.0, 290.0, 475.0, 495.0)],
#               [64, (165.0, 185.0, 175.0, 195.0)],
#               [74, (290.0, 310.0, 130.0, 150.0)],
#               [90, (245.0, 265.0, 160.0, 180.0)],
#               [92, (240.0, 260.0, 310.0, 330.0), (492.0, 512.0, 290.0, 310.0)],
#               [146, (20.0, 40.0, 225.0, 245.0)],
#               [154, (20.0, 40.0, 170.0, 190.0)],
#               [158, (165.0, 185.0, 280.0, 300.0)],
#               [170, (290.0, 310.0, 210.0, 230.0)],
#               [172, (50.0, 70.0, 470.0, 490.0)],
#               [188, (35.0, 55.0, 60.0, 80.0)]]
#
# custom_fig_dir = os.path.join(local_config.fig_dir, '20150816', 'custom_images')
# custom_results_dir = os.path.join(local_config.results_dir, '20150816')
#
#
# def im_coord_given_pos(im_pos):
#     row = im_pos / 20
#     col = im_pos % 20
#     return '%s%d' % (string.uppercase[row], col)
#
# def im_pos_given_idx(im_idx):
#     return int(im_idx / len(nd2.channels)) + 1
#
# def im_coord_given_idx(im_idx):
#     return im_coord_given_pos(im_pos_given_idx(im_idx))
#
# def process(location_idx):
#     location = locations[location_idx]
#     im_idx = location[0]
#     phix_im_idx = im_idx + 1
#     im_coord = im_coord_given_idx(im_idx)
#     all_coords = location[1:]
#     sexcat_fpath = os.path.join(sexcat_dir, str(phix_im_idx) + '.cat')
#
#     pep_im = nd2[im_idx].data
#
#     fic = fastqimagecorrelator.FastqImageAligner(project_name)
#     fic.load_phiX()
#     fic.set_image_data(im=nd2[phix_im_idx].data, objective=60, fpath=im_coord, median_normalize=True)
#     fic.set_sexcat_from_file(sexcat_fpath)
#     fic.align(possible_tile_keys, rotation_est, fq_w_est=fq_w_est, snr_thresh=1.2, min_hits=15)
#
#     all_fic = fastqimagecorrelator.FastqImageAligner(project_name)
#     tile_data=local_config.fastq_tiles_given_read_name_fpath(all_read_name_fpath)
#     all_fic.all_reads_fic_from_aligned_fic(fic, tile_data=tile_data)
#     all_aligned_rcs = all_fic.hitting_tiles[0].aligned_rcs
#
#     fig, ax = plt.subplots(figsize=(15, 15))
#     ax.matshow(pep_im, cmap=plt.get_cmap('Blues'))
#     aligned_rcs = all_fic.hitting_tiles[0].aligned_rcs
#     ax.plot(aligned_rcs[:, 1], aligned_rcs[:, 0], '.', color='darkgoldenrod', alpha=0.6, markersize=3,
#             label='All')
#     aligned_rcs = fic.hitting_tiles[0].aligned_rcs
#     ax.plot(aligned_rcs[:, 1], aligned_rcs[:, 0], 'r.', alpha=0.8, markersize=3, label='PhiX')
#     ax.set_xlim((0, 512))
#     ax.set_ylim((512, 0))
#     ax.set_title('Image %s' % im_coord)
#     fig.savefig(os.path.join(custom_fig_dir, 'im_%s_all_points_overlay.pdf' % (im_coord)))
#
#     def find_points(rlim, clim):
#         out_rcs = []
#         out_read_names = []
#         for i, pt in enumerate(all_aligned_rcs):
#             if rlim[0] <= pt[0] <= rlim[1] and clim[0] <= pt[1] <= clim[1]:
#                 out_rcs.append(pt)
#                 out_read_names.append(all_fic.hitting_tiles[0].read_names[i])
#         return out_rcs, out_read_names
#
#     for i, zoom_coords in enumerate(all_coords):
#         zoom_rcs, zoom_read_names = find_points(zoom_coords[2:], zoom_coords[:2])
#         hot_rcs = []
#         hot_read_names = []
#         for pt, read_name in zip(zoom_rcs, zoom_read_names):
#             r = int(round(pt[0]))
#             start_r = max(0, r-1)
#             c = int(round(pt[1]))
#             start_c = max(0, r-1)
#             if pep_im[start_r:r+2, start_c:c+2].max() > 2:
#                 hot_rcs.append(pt)
#                 hot_read_names.append(read_name)
#         hot_rcs = np.array(hot_rcs)
#
#         fig, ax = plt.subplots(figsize=(15, 15))
#         ax.matshow(pep_im, cmap=plt.get_cmap('Blues'))
#         ax.plot(all_aligned_rcs[:, 1], all_aligned_rcs[:, 0], 'o', color='darkgoldenrod', alpha=0.6,
#                 markersize=10, label='All')
#         ax.plot(aligned_rcs[:, 1], aligned_rcs[:, 0], 'ro', alpha=0.8, markersize=10, label='PhiX')
#         if hot_rcs:
#             ax.plot(hot_rcs[:, 1], hot_rcs[:, 0], 'o', mfc='none', markersize=14, label='Called')
#         ax.set_xlim(zoom_coords[:2])
#         ax.set_ylim(zoom_coords[2:][::-1])
#         ax.set_title('Image %s' % im_coord)
#         ax.legend()
#         fig.savefig(os.path.join(custom_fig_dir, 'im_%s_zoom_point_%d_overlay.pdf' % (im_coord, i)))
#
#         out_fpath = os.path.join(custom_results_dir,
#                                  'im_%s_zoom_point_%d_called_reads.txt' % (im_coord, i))
#         with open(out_fpath, 'w') as out:
#             out.write('\n'.join(hot_read_names))
#
#
#
# if __name__ == '__main__':
#     fmt = '%s <location_idx>' % sys.argv[0]
#     if len(sys.argv) != len(fmt.split()):
#         sys.exit('Usage: ' + fmt)
#     location_idx = int(sys.argv[1])
#     process(location_idx)
