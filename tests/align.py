import unittest
from chimp.align import get_expected_tile_map


class TestAlign(unittest.TestCase):
    def test_get_expected_tile_map_left_edge_nonmaxed_tiles(self):
        tile_map = get_expected_tile_map(4, 14, 0, 119)
        self.assertListEqual(tile_map[0], [4, 5])

    def test_get_expected_tile_map_left_edge_maxed_tiles(self):
        tile_map = get_expected_tile_map(1, 19, 0, 119)
        self.assertListEqual(tile_map[0], [1, 2])

    def test_get_expected_tile_map_right_edge_nonmaxed_tiles(self):
        tile_map = get_expected_tile_map(4, 14, 0, 119)
        self.assertListEqual(tile_map[119], [14, 13])

    def test_get_expected_tile_map_right_edge_maxed_tiles(self):
        tile_map = get_expected_tile_map(1, 19, 0, 119)
        self.assertListEqual(tile_map[119], [19, 18])
