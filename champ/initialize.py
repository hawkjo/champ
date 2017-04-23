import os
import yaml
from champ.error import fail


def save(clargs):
    filename = os.path.join(clargs.image_directory, 'champ.yml')
    with open(filename, 'w+') as f:
        data = {'mapped_reads': os.path.abspath(clargs.mapped_reads),
                'microns_per_pixel': clargs.microns_per_pixel,
                'chip_type': str(clargs.chip),
                'ports_on_right': clargs.ports_on_right,
                'alignment_channel': clargs.alignment_channel,
                'alternate_fiducial_reads': clargs.alternate_fiducial_reads,
                'alternate_perfect_target_reads_filename': clargs.alternate_perfect_target_reads_filename,
                'alternate_good_target_reads_filename': clargs.alternate_good_target_reads_filename,
                'flipud': clargs.flipud,
                'fliplr': clargs.fliplr,
                'perfect_target_name': clargs.perfect_target_name,
                'phix_aligned': False,
                'preprocessed': False,
                'protein_channels_aligned': []}
        yaml.dump(data, f)


def update(image_directory, metadata):
    filename = get_existing_metadata_filename(image_directory)
    with open(filename, 'w+') as f:
        yaml.dump(metadata, f)


def load(image_directory):
    filename = get_existing_metadata_filename(image_directory)
    try:
        with open(filename) as fh:
            return yaml.load(fh)
    except IOError:
        fail("The image directory you provided (%s) has not been initialized. We need you to provide metadata with"
             "the 'champ init' command first." % str(image_directory))
    except:
        fail("Something is wrong with the metadata file in the image directory. Try rerunning 'champ init'.", 2)


def get_existing_metadata_filename(image_directory):
    # Look for the metadata file, and if it doesn't exist, try again with the old extension to allow backwards compatibility
    filename = os.path.join(image_directory, 'champ.yml')
    if not os.path.exists(filename):
        filename = os.path.join(image_directory, 'champ.yaml')
    return filename
