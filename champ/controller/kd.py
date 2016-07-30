import yaml
from champ import intensity


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
    intensity.main(clargs, target_name, target_sequence, off_target_sequence)