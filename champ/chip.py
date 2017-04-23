from collections import defaultdict
import itertools


class BaseChip(object):
    def __init__(self, tile_count):
        self._tile_count = tile_count
        self._lane = 1

    def expected_tile_map(self, left_tiles, right_tiles, min_column, max_column):
        # Creates a dictionary that relates each column of microscope images to its expected tile, +/- 1.
        # Works regardless of whether everything is flipped upsidedown (i.e. the lower tile is on the
        # right side)
        tile_map = defaultdict(list)

        # We gets lists of tiles, so we have to work out the minimum and maximum number in a slightly
        # complicated way
        left_tiles = [int(tile[-4:]) for tile in left_tiles]
        right_tiles = [int(tile[-4:]) for tile in right_tiles]
        min_tile = min(itertools.chain(left_tiles, right_tiles))
        max_tile = max(itertools.chain(left_tiles, right_tiles))

        # Keep track of whether we'll have to invert all the associations of tiles and columns
        invert_map = True if min_tile not in left_tiles else False

        # Find the "tiles per column" factor so we can map a column to a tile
        normalization_factor = float(abs(max_tile - min_tile) + 1) / float(max_column - min_column)

        # Build up the map
        for column in range(min_column, max_column + 1):
            expected = int(round(normalization_factor * column)) - 1
            expected = min(self._tile_count, max(0, expected)) + min_tile
            tile_map_column = column if not invert_map else max_column - column
            # We definitely need to check the predicted tile
            tile_map[tile_map_column].append(self._format_tile_number(expected))
            # If we're at a boundary, we just want to check the adjacent tile towards the middle
            # If we're in the middle, we want to check the tiles on either side
            if expected < max_tile:
                tile_map[tile_map_column].append(self._format_tile_number(expected + 1))
            if expected > min_tile:
                tile_map[tile_map_column].append(self._format_tile_number(expected - 1))
        return tile_map

    def _format_tile_number(self, tile):
        return 'lane{lane}tile{tile}'.format(lane=self._lane, tile=tile)

    @property
    def tile_count(self):
        return self._tile_count


class Miseq(BaseChip):
    def __init__(self, ports_on_right, side=2):
        # ports_on_right means that when imaging, the two fluid inlet ports are towards the right, at least on side 2
        super(Miseq, self).__init__(19)
        self._lower_tiles = [self._format_tile_number(int("%d100" % side) + num) for num in range(1, 11)]
        self._higher_tiles = [self._format_tile_number(int("%d100" % side) + num) for num in reversed(range(11, 20))]
        self._ports_on_right = ports_on_right
        self.tile_width = 935.0
        self.rotation_estimate = 180.0
        self._lane = 1

    def __str__(self):
        return 'miseq'

    @property
    def right_side_tiles(self):
        return self._higher_tiles if not self._ports_on_right else self._lower_tiles

    @property
    def left_side_tiles(self):
        return self._lower_tiles if not self._ports_on_right else self._higher_tiles


class Hiseq(BaseChip):
    def __init__(self, ports_on_right, lane=4, side=2):
        super(Hiseq, self).__init__(28)
        self._lower_tiles = [self._format_tile_number(int("%d100" % side) + num) for num in range(1, 14)]
        self._higher_tiles = [self._format_tile_number(int("%d100" % side) + num) for num in reversed(range(15, 28))]
        self._ports_on_right = ports_on_right
        self._lane = lane
        self.tile_width = 1097.0
        self.rotation_estimate = 180.0

    def __str__(self):
        return 'hiseq'

    @property
    def right_side_tiles(self):
        return self._higher_tiles if not self._ports_on_right else self._lower_tiles

    @property
    def left_side_tiles(self):
        return self._lower_tiles if not self._ports_on_right else self._higher_tiles


def load(chip_name):
    chips = {'miseq': Miseq,
             'hiseq': Hiseq}
    return chips[chip_name]
