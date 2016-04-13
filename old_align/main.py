from fastqimagealigner import FastqImageAligner
import os
from chimp import sextraction


def load_sextraction(experiment, nd2_name, image_index):
    with open(experiment.get_sexcat_path(nd2_name, image_index)) as f:
        return sextraction.Sextraction(f)


def align(image, tile, sexcat, params):
    """

    :param image:
    :param tile:
    :param params:
    :return:
    """
    fic = FastqImageAligner(params.chip_id)
    # fic.set_image_data(im=nd2[image_index], objective=params.objective,
    #                    fpath=str(image_index), median_normalize=True)
    # fic.sexcat = load_sextraction(experiment, nd2_name, image_index)

    fic.align(possible_tile_keys,
              alignment_parameters.rotation_estimate,
              alignment_parameters.fq_w_est,
              snr_thresh=alignment_parameters.snr_threshold,
              min_hits=alignment_parameters.min_hits,
              hit_type=('exclusive', 'good_mutual'))
