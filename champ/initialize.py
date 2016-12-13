import os
import yaml
from champ.error import fail


def save(clargs):
    metadata_filename = os.path.join(clargs.image_directory, 'champ.yml')
    cache_filename = os.path.join(clargs.image_directory, 'cache.yml')

    # save static data that should never change unless the user makes a mistake
    with open(metadata_filename, 'w') as f:
        data = {'chip_name': clargs.chip_name,
                'mapped_reads': os.path.abspath(clargs.mapped_reads),
                'parsed_reads': os.path.abspath(clargs.parsed_reads),
                'lda_weights': os.path.abspath(clargs.nonneg_lda_weights_path),
                'microns_per_pixel': clargs.microns_per_pixel,
                'chip_type': str(clargs.chip),
                'ports_on_right': clargs.ports_on_right,
                'alignment_channel': clargs.alignment_channel,
                'alternate_fiducial_reads': clargs.alternate_fiducial_reads,
                'alternate_perfect_target_reads_filename': clargs.alternate_perfect_target_reads_filename,
                'alternate_good_target_reads_filename': clargs.alternate_good_target_reads_filename,
                'flipud': clargs.flipud,
                'fliplr': clargs.fliplr,
                'perfect_target_name': clargs.perfect_target_name
                }
        yaml.dump(data, f)

    # save empty cache
    with open(cache_filename, 'w') as f:
        data = {'phix_aligned': False,
                'preprocessed': False,
                'protein_channels_aligned': []}
        yaml.dump(data, f)


def update_cache(image_directory, cache):
    filename = os.path.join(image_directory, 'cache.yml')
    with open(filename, 'w') as f:
        yaml.dump(cache, f)


def load_cache_and_metadata(image_directory):
    metadata_filename = get_existing_metadata_filename(image_directory)
    cache_filename = os.path.join(image_directory, 'cache.yml')
    try:
        with open(metadata_filename) as fh:
            metadata = yaml.load(fh)
    except IOError:
        fail("The image directory you provided (%s) has not been initialized. We need you to provide metadata with"
             "the 'champ init' command first." % str(image_directory))
    except:
        fail("Something is wrong with the metadata file in the image directory. Try rerunning 'champ init'.", 2)
    else:
        try:
            with open(cache_filename) as fh:
                cache = yaml.load(fh)
        except IOError:
            # make separation of cache file backwards compatible with datasets that didn't start with it
            cache = {}
            cache['phix_aligned'] = metadata['phix_aligned']
            cache['preprocessed'] = metadata['preprocessed']
            cache['protein_channels_aligned'] = metadata['protein_channels_aligned']
            if 'end_tiles' in metadata:
                cache['end_tiles'] = metadata['end_tiles']
    return metadata, cache


def get_existing_metadata_filename(image_directory):
    # Look for the metadata file, and if it doesn't exist,
    # try again with the old extension to allow backwards compatibility
    filename = os.path.join(image_directory, 'champ.yml')
    if not os.path.exists(filename):
        filename = os.path.join(image_directory, 'champ.yaml')
    return filename
