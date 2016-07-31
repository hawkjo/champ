import yaml
from champ import intensity, initialize


class TargetInfo(object):
    """
    This thing reads the target info file and provides methods to access various targets
    """
    def __init__(self, path, on_target_label, off_target_label):
        self.on_target_label = on_target_label.lower()
        self.off_target_label = off_target_label.lower()
        with open(path) as f:
            self._data = {key.lower(): val for key, val in yaml.load(f).items()}

    def on_target_sequence(self):
        return self._data[self.on_target_label]

    def off_target_sequence(self):
        return self._data[self.off_target_label]


def main(clargs):
    metadata = initialize.load(clargs.image_directory)
    target_info = TargetInfo(clargs.target_data_file, clargs.target_label, clargs.off_target_label)
    intensity.main(metadata, clargs.image_directory, target_info)