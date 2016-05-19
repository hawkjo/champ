from chimp.align import get_expected_tile_map
from pprint import pprint


class Tile(object):
    def __init__(self, key):
        self.key = key


pprint(dict(get_expected_tile_map([Tile('2104')], [Tile('2114')], 0, 58)))