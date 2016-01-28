import matplotlib
matplotlib.use('agg')
import sys
import os
import fastqimagealigner
import local_config
import nd2reader
import nd2tools
import logging
import reads
import params

log = logging.getLogger(__name__)


def fast_possible_tile_keys(nd2, im_idx, min_tile, max_tile):
    coord_info, xs, ys, zs, pos_names, rows, cols = nd2tools.get_nd2_image_coord_info(nd2)
    cols.sort()
    pos_name = nd2tools.convert_nd2_coordinates(nd2, im_idx=im_idx, outfmt='pos_name')
    col_idx = cols.index(pos_name[1:])
    expected_tile = int(min_tile + col_idx * float(max_tile - min_tile)/(len(cols)-1))
    return tile_keys_given_nums(range(expected_tile - 1, expected_tile + 2))


def tile_keys_given_nums(tile_nums):
    return ['lane1tile{0}'.format(tile_num) for tile_num in tile_nums]


def process_fig(align_run_name, base_directory, nd2_fpath, align_param_fpath, im_idx):
    file_structure = local_config.FileStructure(base_directory)
    im_idx = int(im_idx)
    alignment_parameters = params.get_align_params(align_param_fpath)
    nd2 = nd2reader.Nd2(nd2_fpath)
    if alignment_parameters.strategy == 'fast':
        possible_tile_keys = fast_possible_tile_keys(nd2, im_idx, alignment_parameters.min_tile_num, alignment_parameters.max_tile_num)
    else:
        assert alignment_parameters.strategy == 'slow'
        possible_tile_keys = tile_keys_given_nums(range(alignment_parameters.min_tile_num, alignment_parameters.max_tile_num + 1))
    bname = os.path.splitext(os.path.basename(nd2_fpath))[0]
    sexcat_fpath = os.path.join(os.path.splitext(nd2_fpath)[0], '%d.cat' % im_idx)
    
    fig_dir = os.path.join(file_structure.figure_directory, align_run_name, bname)
    results_dir = os.path.join(base_directory, 'results', align_run_name, bname)
    for d in (fig_dir, results_dir):
        if not os.path.exists(d):
            os.makedirs(d)
    all_read_rcs_fpath = os.path.join(results_dir, '{}_all_read_rcs.txt'.format(im_idx))

    if os.path.isfile(all_read_rcs_fpath):
        log.debug('%s, %s already done.' % (bname, im_idx))

    intensity_fpath = os.path.join(results_dir, '{}_intensities.txt'.format(im_idx))
    stats_fpath = os.path.join(results_dir, '{}_stats.txt'.format(im_idx))
    fic = fastqimagealigner.FastqImageAligner(alignment_parameters.project_name, file_structure)
    tile_data = reads.get_read_names(alignment_parameters.aligning_read_names_fpath)
    fic.load_reads(tile_data)
    fic.set_image_data(im=nd2[im_idx], objective=alignment_parameters.objective, fpath=str(im_idx), median_normalize=True)
    fic.set_sexcat_from_file(sexcat_fpath)
    fic.align(possible_tile_keys, alignment_parameters.rotation_est, alignment_parameters.fq_w_est,
              snr_thresh=alignment_parameters.snr_threshold,
              min_hits=alignment_parameters.min_hits,
              hit_type=['exclusive', 'good_mutual'])
    log.debug("%s %s %s %s" % (alignment_parameters.project_name, bname, im_idx, ','.join(tile.key for tile in fic.hitting_tiles)))
    
    fic.output_intensity_results(intensity_fpath)
    fic.write_alignment_stats(stats_fpath)

    ax = fic.plot_all_hits()
    ax.figure.savefig(os.path.join(fig_dir, '{}_all_hits.pdf'.format(im_idx)))

    ax = fic.plot_hit_hists()
    ax.figure.savefig(os.path.join(fig_dir, '{}_hit_hists.pdf'.format(im_idx)))

    all_fic = fastqimagealigner.FastqImageAligner(alignment_parameters.project_name, file_structure)
    tile_data = reads.get_read_names(alignment_parameters.all_read_names_fpath)
    all_fic.all_reads_fic_from_aligned_fic(fic, tile_data)
    all_fic.write_read_names_rcs(all_read_rcs_fpath)

if __name__ == '__main__':
    fmt = '{0} <align_run_name> <nd2_fpath> <align_param_file> <im_idx>'.format(sys.argv[0])
    if len(sys.argv) != len(fmt.split()):
        sys.exit('Usage: ' + fmt)
    process_fig(*sys.argv[1:])
