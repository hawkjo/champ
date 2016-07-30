import yaml
from champ import intensity, initialize


class TargetInfo(object):
    """
    This thing reads the target info file and provides methods to access various targets
    """
    def __init__(self, path, on_target, off_target):
        pass

    def on_target(self):
        pass

    def off_target(self):
        pass


def main(clargs):
    metadata = initialize.load(clargs.image_directory)
    target_info = TargetInfo(clargs.target_data_file, clargs.target_label, clargs.off_target_label)
    intensity.main(metadata, clargs.image_directory, target_info)