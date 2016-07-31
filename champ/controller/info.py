import logging
from champ import projectinfo


log = logging.getLogger(__name__)


def main(clargs):
    channels = projectinfo.load_channels(clargs.image_directory)
    print("Channels:")
    for channel in channels:
        print(channel)
