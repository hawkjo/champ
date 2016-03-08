from chimp.model.sextractor import Sextraction
from chimp.model.grid import GridImages
from chimp.model import tile
from chimp import fastq
from functools import partial
from nd2reader import Nd2
import numpy as np
import os


def load_sexcat(directory, index):
    with open(os.path.join(directory, '%s.cat' % index)) as f:
        return Sextraction(f)


def find_end_tiles():
    """
    What's the deliverable here?

    two numerical tile numbers and two column numbers

    so we just left_iter and right_iter until we find a good alignment
    now we know two tiles and we know the images that fall inside of them


    """
    pass


def pad_images(tile_image, microscope_image):
    """
    Pad to 4096 pixels, or whatever.

    """
    pad = np.abs(tile_image.shape - microscope_image.shape)
    return pad_image(tile_image, pad), pad_image(microscope_image, pad)


def pad_image(image, pad):
    return np.pad(image, ((0, pad[0]), (0, pad[1])), mode='constant')


def main(arguments, image_files):
    """
    Let's think through this.
    So we have like 11 different concentrations.
    The tiles are always the same (data-wise) though each concentration might have reads in tiles that other's don't.

    The microscope image and sexcat are different for each grid in each concentration.

    So we have 19 tiles, and 7*60*11 grid images.

    What's the output?

    Alignment (transform, rotate, scale) for EVERY image in EVERY concentration, where alignment occurred

    So we want an AlignmentResults object for each concentration. We do that so we can checkpoint at the end of each,
    because this will fail. Need to build restart ability eventually.

    END FINDING
    Take each image using left and right iter
    Try to align to 1-4 or 19-16 and also controls
    Get the outermost tiles and the columns of the images

    Do simple linear interpolation for tiles/columns
    Do precision alignment for every image inbetween the bounding columns, as you iterate over interpolated tiles
    Put results into the AlignmentResults object

    CONSIDER
    For rough alignment, you need microscope image and tile
    For precision alignment, you need sextractor image and tile

    Not the output:
    Intensity data (from sextractor)
    Read names with real coordinates

    """
    read_data = fastq.load_mapped_reads(arguments.alignment_reads,
                                        arguments.ignore_side_1)
    tiles = tile.load_tile_manager(nd2.pixel_microns, read_data)



# def main(base_image_name, alignment_channel=None, alignment_offset=None):
#     nd2 = Nd2('%s.nd2' % base_image_name)
#     loader = partial(load_sexcat, base_image_name)
#     grid = GridImages(nd2, loader, alignment_channel, alignment_offset)
#     md = grid.get(5, 5)
#     io.imsave('jim.pdf', md.seimage)


# if __name__ == '__main__':
#     main('/var/experiments/151118/15-11-18_SA15243_Cascade-TA_1nM-007', alignment_offset=1)
