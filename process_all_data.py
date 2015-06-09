import sys
import os
import fastqimagecorrelator
import matplotlib.pyplot as plt
from pathos.multiprocessing import ProcessingPool
import local_config


def get_align_params(align_param_fpath):
    #d = {name: value for line in open(align_param_fpath) for name, value in line.strip().split()}
    d = {}
    for line in open(align_param_fpath):
        name, value = line.strip().split()
        d[name] = value

    possible_tile_keys = ['lane1tile{0}'.format(tile_num) 
                          for tile_num in map(int, d['possible_tiles'].split(','))]

    return (d['project_name'],
            int(d['objective']),
            float(d['rotation_est']),
            float(d['fq_w_est']),
            possible_tile_keys,
           )


def process_fig(align_run_name, align_param_fpath, im_fpath):
    project_name, objective, rotation_est, fq_w_est, possible_tile_keys \
            = get_align_params(align_param_fpath)
    bname = os.path.splitext(os.path.basename(im_fpath))[0]
    sexcat_fpath = im_fpath.replace('.tif', '.cat')
    
    fic = fastqimagecorrelator.FastqImageCorrelator(project_name)
    fic.load_phiX()
    fic.set_image_data(im_fpath, objective, median_normalize=True)
    fic.set_sexcat_from_file(sexcat_fpath)
    fic.align(possible_tile_keys, rotation_est, fq_w_est)
    print project_name, bname, ','.join(tile.key for tile in fic.hitting_tiles)
    
    fig_dir = os.path.join(local_config.fig_dir, align_run_name)
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    ax = fic.plot_all_hits()
    ax.figure.savefig(os.path.join(fig_dir, '{0}_all_hits.pdf'.format(bname)))

    ax = fic.plot_hit_hists()
    ax.figure.savefig(os.path.join(fig_dir, '{0}_hit_hists.pdf'.format(bname)))

    results_dir = os.path.join(local_config.base_dir, 'results', align_run_name)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    intensity_fpath = os.path.join(results_dir, '{0}_intensities.txt'.format(bname))
    stats_fpath = os.path.join(results_dir, '{0}_stats.txt'.format(bname))
    fic.output_intensity_results(intensity_fpath)
    fic.write_alignment_stats(stats_fpath)

    all_fic = fastqimagecorrelator.FastqImageCorrelator(project_name)
    all_fic.all_reads_fic_from_aligned_fic(fic)
    all_read_rcs_fpath = os.path.join(results_dir, '{0}_all_read_rcs.txt'.format(bname))
    all_fic.write_read_names_rcs(all_read_rcs_fpath)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit('Usage: {0} <align_run_name> <align_param_file> <image_file>'.format(sys.argv[0]))
    process_fig(*sys.argv[1:])
