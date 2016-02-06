import sys
import os
import fastqimagealigner
import config
import nd2reader
import nd2tools
import logging
import reads

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


def process_fig(base_directory, chip_id, strategy, nd2_filename, im_idx):
    assert strategy in ('fast', 'slow'), 'Invalid alignment strategy'
    file_structure = config.Experiment(base_directory)
    im_idx = int(im_idx)
    alignment_parameters = config.AlignmentParameters(base_directory, chip_id)
    nd2 = nd2reader.Nd2('{base_directory}{sep}{nd2_filename}'.format(base_directory=base_directory,
                                                                     nd2_filename=nd2_filename,
                                                                     sep=os.path.sep))
    if strategy == 'fast':
        # +4 and -5 are temporary hacks to see if something works
        # if you are on the master branch and you're reading this, Jim is a failure as a programmer, and as a human being
        # remind him that he has dishonored not only himself, but his children, and his children's children, for nine generations to come
        possible_tile_keys = fast_possible_tile_keys(nd2, im_idx, alignment_parameters.min_tile_num + 4, alignment_parameters.max_tile_num - 5)
    else:
        possible_tile_keys = tile_keys_given_nums(range(alignment_parameters.min_tile_num, alignment_parameters.max_tile_num + 1))

    for directory in (file_structure.figure_directory, file_structure.results_directory):
        if not os.path.exists(directory):
            os.makedirs(directory)
    all_read_rcs_filepath = os.path.join(file_structure.results_directory, '{}_all_read_rcs.txt'.format(im_idx))
    base_nd2_name = os.path.splitext(os.path.basename(nd2_filename))[0]
    log.debug("bname %s" % base_nd2_name)
    if os.path.isfile(all_read_rcs_filepath):
        log.debug('%s, %s already done.' % (base_nd2_name, im_idx))

    sexcat_fpath = os.path.join(base_directory, os.path.splitext(nd2_filename)[0], '%d.cat' % im_idx)
    intensity_fpath = os.path.join(file_structure.results_directory, '{}_intensities.txt'.format(im_idx))
    stats_fpath = os.path.join(file_structure.results_directory, '{}_stats.txt'.format(im_idx))
    fic = fastqimagealigner.FastqImageAligner(chip_id, file_structure)
    tile_data = reads.get_read_names(alignment_parameters.aligning_read_names_filepath)
    fic.load_reads(tile_data)
    fic.set_image_data(im=nd2[im_idx], objective=alignment_parameters.objective, fpath=str(im_idx), median_normalize=True)
    fic.set_sexcat_from_file(sexcat_fpath)
    fic.align(possible_tile_keys,
              alignment_parameters.rotation_estimate,
              alignment_parameters.fq_w_est,
              snr_thresh=alignment_parameters.snr_threshold,
              min_hits=alignment_parameters.min_hits,
              hit_type=('exclusive', 'good_mutual'))
    log.debug("%s %s %s %s" % (alignment_parameters.chip_id, base_nd2_name, im_idx, ','.join(tile.key for tile in fic.hitting_tiles)))
    
    fic.output_intensity_results(intensity_fpath)
    fic.write_alignment_stats(stats_fpath)

    ax = fic.plot_all_hits()
    ax.figure.savefig(os.path.join(file_structure.figure_directory, '{}_all_hits.pdf'.format(im_idx)))

    ax = fic.plot_hit_hists()
    ax.figure.savefig(os.path.join(file_structure.figure_directory, '{}_hit_hists.pdf'.format(im_idx)))

    all_fic = fastqimagealigner.FastqImageAligner(alignment_parameters.chip_id, file_structure)
    tile_data = reads.get_read_names(alignment_parameters.all_read_names_filepath)
    all_fic.all_reads_fic_from_aligned_fic(fic, tile_data)
    all_fic.write_read_names_rcs(all_read_rcs_filepath)

if __name__ == '__main__':
    fmt = '{0} <base_directory> <chip_id> <strategy> <nd2_fpath> <im_idx>'.format(sys.argv[0])
    if len(sys.argv) != len(fmt.split()):
        sys.exit('Usage: ' + fmt)
    process_fig(*sys.argv[1:])
