import logging
from champ import initialize


log = logging.getLogger(__name__)


def main(clargs):
    """ Stores details of the experiment that are needed for all analyses. This accomplishes two things: first, it
        reduces the number of arguments that have to be specified in further commands, and also acts as a way of
        documenting the experiment. """
    channels = initialize.determine_channel_names(clargs.image_directory)
    if len(channels) > 1 and not clargs.alignment_channel:
        alignment_channel = initialize.request_alignment_channel(channels)
    initialize.save_metadata(clargs)
