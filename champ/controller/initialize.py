import logging
from champ import initialize


log = logging.getLogger(__name__)


def main(clargs):
    """ Stores details of the experiment that are needed for all analyses. This accomplishes two things: first, it
        reduces the number of arguments that have to be specified in further commands, and also acts as a way of
        documenting the experiment. """
    initialize.save_metadata(clargs)
