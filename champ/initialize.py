import yaml
import os
from error import fail


def save(clargs):
    filename = os.path.join(clargs.image_directory, 'champ.yaml')
    with open(filename, 'w+') as f:
        data = {'chip_name': clargs.chip_name,
                'mapped_reads': os.path.abspath(clargs.mapped_reads),
                'parsed_reads': os.path.abspath(clargs.parsed_reads),
                'lda_weights': os.path.abspath(clargs.nonneg_lda_weights_path),
                'microns_per_pixel': clargs.microns_per_pixel,
                'chip_type': str(clargs.chip),
                'ports_on_right': clargs.ports_on_right,
                'alignment_channel': clargs.alignment_channel,
                'flipud': clargs.flipud,
                'fliplr': clargs.fliplr}
        yaml.dump(data, f)


def load(image_directory):
    filename = os.path.join(image_directory, 'champ.yaml')
    try:
        with open(filename) as fh:
            return yaml.load(fh)
    except IOError:
        fail("The image directory you provided (%s) has not been initialized. We need you to provide metadata with"
             "the 'champ init' command first." % str(image_directory))
    except:
        fail("Something is wrong with the metadata file in the image directory. Try rerunning 'champ init'.", 2)
