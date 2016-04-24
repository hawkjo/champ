import fastqimagealigner
import logging
import os

log = logging.getLogger(__name__)


def process_fig(alignment_parameters, base_name, tile_data, um_per_pixel, experiment, image, possible_tile_keys):
    for directory in (experiment.figure_directory, experiment.results_directory):
        full_directory = os.path.join(directory, base_name)
        if not os.path.exists(full_directory):
            os.makedirs(full_directory)

    sexcat_fpath = os.path.join(base_name, '%d.cat' % image.index)
    fic = fastqimagealigner.FastqImageAligner(experiment)
    fic.load_reads(tile_data)
    fic.set_image_data(image, um_per_pixel)
    fic.set_sexcat_from_file(sexcat_fpath)
    fic.rough_align(possible_tile_keys,
                    alignment_parameters.rotation_estimate,
                    alignment_parameters.fq_w_est,
                    snr_thresh=alignment_parameters.snr_threshold)
    return fic


def write_output(image_index, base_name, fastq_image_aligner, experiment, tile_data):
    intensity_filepath = os.path.join(experiment.results_directory, base_name, '{}_intensities.txt'.format(image_index))
    stats_filepath = os.path.join(experiment.results_directory, base_name, '{}_stats.txt'.format(image_index))
    all_read_rcs_filepath = os.path.join(experiment.results_directory, base_name, '{}_all_read_rcs.txt'.format(image_index))

    fastq_image_aligner.output_intensity_results(intensity_filepath)
    fastq_image_aligner.write_alignment_stats(stats_filepath)
    all_fastq_image_aligner = fastqimagealigner.FastqImageAligner(experiment)
    all_fastq_image_aligner.all_reads_fic_from_aligned_fic(fastq_image_aligner, tile_data)
    all_fastq_image_aligner.write_read_names_rcs(all_read_rcs_filepath)
