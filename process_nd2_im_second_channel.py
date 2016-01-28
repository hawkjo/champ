import matplotlib
matplotlib.use('agg')
import sys
import os
import fastqimagealigner
import local_config
import nd2reader
import logging
import reads

log = logging.getLogger(__name__)


class AlignmentParameters(object):
    def __init__(self, lines):
        self._data = {}
        for line in map(str.strip, lines):
            if not line:
                continue
            name, value = line.split()
            self._data[name] = value

        self._min_hits = self._data.get('min_hits', 15)
        self._strategy = self._data.get('strategy', 'slow')
        assert self._strategy in ['fast', 'slow'], 'Invalid alignment strategy: {strategy}'.format(strategy=self._strategy)

    @property
    def project_name(self):
        return self._data['project_name']

    @property
    def aligning_read_names_fpath(self):
        return self._data['aligning_read_names_fpath']

    @property
    def all_read_names_fpath(self):
        return self._data['all_read_names_fpath']

    @property
    def objective(self):
        return int(self._data['objective'])

    @property
    def aligned_im_idx_offset(self):
        return int(self._data['aligned_im_idx_offset'])

    @property
    def min_hits(self):
        return int(self._min_hits)


def get_align_params(align_param_fpath):
    with open(align_param_fpath) as lines:
        return AlignmentParameters(lines)


def process_fig(align_run_name, base_directory, nd2_fpath, align_param_fpath, im_idx):
    file_structure = local_config.FileStructure(base_directory)
    im_idx = int(im_idx)
    alignment_parameters = get_align_params(align_param_fpath)
    nd2 = nd2reader.Nd2(nd2_fpath)
    bname = os.path.splitext(os.path.basename(nd2_fpath))[0]
    aligned_im_idx = im_idx + alignment_parameters.aligned_im_idx_offset
    sexcat_fpath = os.path.join(os.path.splitext(nd2_fpath)[0], '%d.cat' % im_idx)

    fig_dir = os.path.join(file_structure.figure_directory, align_run_name, bname)
    results_dir = os.path.join(base_directory, 'results', align_run_name, bname)
    for d in [fig_dir, results_dir]:
        if not os.path.exists(d):
            os.makedirs(d)

    fic = fastqimagealigner.FastqImageAligner(alignment_parameters.project_name, file_structure)
    tile_data = reads.get_read_names(alignment_parameters.aligning_read_names_fpath)
    fic.load_reads(tile_data)
    fic.set_image_data(im=nd2[im_idx], objective=alignment_parameters.objective, fpath=str(im_idx), median_normalize=True)
    fic.set_sexcat_from_file(sexcat_fpath)

    stats_fpath = os.path.join(results_dir, '{}_stats.txt'.format(aligned_im_idx))
    fic.alignment_from_alignment_file(stats_fpath)
    fic.precision_align_only(min_hits=alignment_parameters.min_hits)
    log.debug("%s %s %s %s" % (alignment_parameters.project_name, bname, im_idx, ','.join(tile.key for tile in fic.hitting_tiles)))
    
    intensity_fpath = os.path.join(results_dir, '{}_intensities.txt'.format(im_idx))
    stats_fpath = os.path.join(results_dir, '{}_stats.txt'.format(im_idx))
    fic.output_intensity_results(intensity_fpath)
    fic.write_alignment_stats(stats_fpath)

    ax = fic.plot_all_hits()
    ax.figure.savefig(os.path.join(fig_dir, '{}_all_hits.pdf'.format(im_idx)))

    ax = fic.plot_hit_hists()
    ax.figure.savefig(os.path.join(fig_dir, '{}_hit_hists.pdf'.format(im_idx)))

    all_fic = fastqimagealigner.FastqImageAligner(alignment_parameters.project_name, file_structure)
    tile_data = reads.get_read_names(alignment_parameters.all_read_names_fpath)
    all_fic.all_reads_fic_from_aligned_fic(fic, tile_data)
    all_read_rcs_fpath = os.path.join(results_dir, '{}_all_read_rcs.txt'.format(im_idx))
    all_fic.write_read_names_rcs(all_read_rcs_fpath)

if __name__ == '__main__':
    fmt = '{0} <align_run_name> <nd2_fpath> <align_param_file> <im_idx>'.format(sys.argv[0])
    if len(sys.argv) != len(fmt.split()):
        sys.exit('Usage: ' + fmt)
    process_fig(*sys.argv[1:])
