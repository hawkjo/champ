from chimp.model.sextractor import Sextraction
from chimp import files


def load_sexcat(path):
    with open(path) as f:
        sextraction = Sextraction(f)
