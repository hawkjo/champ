# Load the read_name_intensities from disk
# Know the delta_y of saturated cluster (passed into function)
# Get pileups
# make the frozenset of read_names
# see if that frozenset is in the LRUCache, which has maxsize=20

# if not, find the kd like normal with fit_one_group_kd
# add it to the LRUCache

# if so, just look up the kd in the kds LRUCache
# add that kd to the position_kds dictionary (which is position_kds[contig][base])

# save positions_kds to disk
# this replaces the first two steps in main_gaff
