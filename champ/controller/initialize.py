import logging
from champ import initialize, error


log = logging.getLogger(__name__)


def main(clargs):
    """ Stores details of the experiment that are needed for all analyses. This accomplishes two things: first, it
        reduces the number of arguments that have to be specified in further commands, and also acts as a way of
        documenting the experiment. """
    channels = initialize.determine_channel_names(clargs.image_directory)
    alignment_channel = clargs.alignment_channel
    if len(channels) > 1 and not clargs.alignment_channel:
        alignment_channel = initialize.request_alignment_channel(channels)
    elif clargs.alignment_channel and clargs.alignment_channel not in channels:
        error.fail("The given alignment channel ('%s') does not exist in the image data. Available channels: %s" % (clargs.alignment_channel, ", ".join(channels)))
    log.debug("Initializing with alignment channel: %s" % alignment_channel)
    initialize.save_metadata(clargs, alignment_channel)
