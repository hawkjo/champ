import logging
import yaml
import os


log = logging.getLogger(__name__)


def main(clargs):
    with open(os.path.join(clargs.image_directory, 'champ.yaml'), 'w+') as f:
        data = {'chip_name': clargs.chip_name,
                'mapped_reads': os.path.abspath(clargs.mapped_reads_dir),
                'microns_per_pixel': clargs.microns_per_pixel,
                'chip_type': clargs.chip_type,
                'ports_on_right': clargs.ports_on_right,
                'alignment_channel': clargs.alignment_channel}
        yaml.dump(data, f)
