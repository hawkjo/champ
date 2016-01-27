import matplotlib
matplotlib.use('agg')
import sys
import os
import fastqimagecorrelator
import matplotlib.pyplot as plt
from pathos.multiprocessing import ProcessingPool
import local_config
import tifffile
import misc


def process_fig(align_run_name, im_fpath):
    im_fname = os.path.basename(im_fpath)
    bname = os.path.splitext(im_fname)[0]
    aligning_read_names_fpath = im_fpath.replace(im_fname, '../read_names/' + bname + '.txt')
    all_read_names_fpath = aligning_read_names_fpath
    sexcat_fpath = im_fpath.replace('.tif', '.cat')

    project_name = 'SA15161'
    objective = 60
    rotation_est = 179.7
    fq_w_est = 936.0
    possible_tile_keys = ['lane1tile2108']
    
    fic = fastqimagecorrelator.FastqImageAligner(project_name)
    tile_data=local_config.fastq_tiles_given_read_name_fpath(aligning_read_names_fpath)
    fic.load_reads(tile_data)
    im = tifffile.imread(im_fpath)
    fic.set_image_data(im=im, fpath=im_fpath, objective=objective)
    fic.set_sexcat_from_file(sexcat_fpath)
    fic.align(possible_tile_keys, rotation_est, fq_w_est, hit_type=['exclusive', 'good_mutual'])
    print project_name, bname, ','.join(tile.key for tile in fic.hitting_tiles)
    
    fig_dir = os.path.join(local_config.fig_dir, align_run_name)
    results_dir = os.path.join(local_config.base_dir, 'results', align_run_name)
    for d in [fig_dir, results_dir]:
        if not os.path.exists(d):
            os.makedirs(d)

    ax = fic.plot_all_hits()
    ax.figure.savefig(os.path.join(fig_dir, '{}_all_hits.pdf'.format(bname)))

    ax = fic.plot_hit_hists()
    ax.figure.savefig(os.path.join(fig_dir, '{}_hit_hists.pdf'.format(bname)))

    intensity_fpath = os.path.join(results_dir, '{}_intensities.txt'.format(bname))
    stats_fpath = os.path.join(results_dir, '{}_stats.txt'.format(bname))
    fic.output_intensity_results(intensity_fpath)
    fic.write_alignment_stats(stats_fpath)

    all_fic = fastqimagecorrelator.FastqImageAligner(project_name)
    tile_data = local_config.fastq_tiles_given_read_name_fpath(all_read_names_fpath)
    all_fic.all_reads_fic_from_aligned_fic(fic, tile_data)
    all_read_rcs_fpath = os.path.join(results_dir, '{}_all_read_rcs.txt'.format(bname))
    all_fic.write_read_names_rcs(all_read_rcs_fpath)

if __name__ == '__main__':
    fmt = '{0} <align_run_name> <im_fpath>'.format(sys.argv[0])
    if len(sys.argv) != len(fmt.split()):
        sys.exit('Usage: ' + fmt)
    process_fig(*sys.argv[1:])
