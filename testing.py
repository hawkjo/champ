import sys
import os
import random
import time
import scipy.io
import scipy.spatial
import scipy.stats
from sklearn import mixture
from sklearn import neighbors
import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
from local_config import base_dir

import matplotlib
sys.path.append('/home/jah/projects/sequencing/code/')
from Sequencing import Visualize

#--------------------------------------------------------------------------------
# Run params
#--------------------------------------------------------------------------------

# The fraction of the subsample not to delete
frac_samp_seen = 0.95

# z = sigma / (2 sqrt(rho))
# sigma*eye(2) is the covariance matrix of the measurement noise
# rho is the expected points per unit area of the poisson process
# Hence, 1/(2 sqrt(rho)) is the expected distance to nearest neighbor.
z = 0.1

# Mapping variables
# s is the x-parity: 1=no flip, -1=lrflip
s = -1
scale = 10
theta = 2
offset = np.array([10000, 20000])

#--------------------------------------------------------------------------------
# Read in tiles and find basic size
#--------------------------------------------------------------------------------

tiles = scipy.io.loadmat(os.path.join(base_dir, 'example_data', 'locs.mat'))

ref_points = tiles.values()[0]

x_min = 1000
x_max = 1000
y_min = 1000
y_max = 1000

for name, tile in tiles.items():
    if 'tile' not in name:
        continue
    print name
    tmp_x_min = min(tile[:,0])
    tmp_x_max = max(tile[:,0])
    tmp_y_min = min(tile[:,1])
    tmp_y_max = max(tile[:,1])
    print 'x extrema:', tmp_x_min, tmp_x_max
    print 'y extrema:', tmp_y_min, tmp_y_max
    print tile.shape

    if tmp_x_min < x_min:
        x_min = tmp_x_min
    if tmp_x_max > x_max:
        x_max = tmp_x_max
    if tmp_y_min < y_min:
        y_min = tmp_y_min
    if tmp_y_max > y_max:
        y_max = tmp_y_max
#
#    hull = scipy.spatial.ConvexHull(tile)
#    plt.plot(tile[hull.vertices,0], tile[hull.vertices,1])
#plt.show()
x_min = int(round(x_min))
x_max = int(round(x_max))
y_min = int(round(y_min))
y_max = int(round(y_max))

print
print 'x range: %d-%d' % (x_min, x_max)
print 'y range: %d-%d' % (y_min, y_max)

x_mid = float(x_min + x_max)/2
y_mid = float(y_min + y_max)/2
center = np.array([x_mid, y_mid])
r = float(x_max - x_min)/4
print 'center:', center
print 'r:', r

#--------------------------------------------------------------------------------
# Build k-d tree and get raw sample of points
#--------------------------------------------------------------------------------
ref_tree = scipy.spatial.KDTree(ref_points)
ref_samp_points = ref_points[ref_tree.query_ball_point(center, r)]

print
print 'ref_samp_points:'
print ref_samp_points
print ref_samp_points.shape

#--------------------------------------------------------------------------------
# Delete some points
#--------------------------------------------------------------------------------

num_to_keep = int(frac_samp_seen*ref_samp_points.shape[0])
keepers = np.arange(ref_samp_points.shape[0])
ref_samp_points = ref_samp_points[np.random.choice(keepers, size=num_to_keep, replace=False)]
print
print 'num_to_keep:', num_to_keep
print 'ref_samp_points:'
print ref_samp_points
print ref_samp_points.shape

#--------------------------------------------------------------------------------
# Add noise to points
#--------------------------------------------------------------------------------

rho = float(ref_samp_points.shape[0])/(np.pi*r**2)
sigma = z * 2 * np.sqrt(rho)
mean = np.zeros((2,))
cov = sigma**2 * np.eye(2)
samp_points = ref_samp_points + np.random.multivariate_normal(mean, cov, num_to_keep)

print
print 'rho:', rho
print '1/(2 sqrt(rho)):', 1.0/(2*np.sqrt(rho))
print 'z:', z
print 'sigma:', sigma
print 'mean:', mean
print 'cov:'
print cov
print 'samp_points:'
print samp_points
print samp_points.shape

#--------------------------------------------------------------------------------
# Scale, rotate, translate, flip
#--------------------------------------------------------------------------------

def right_2d_rotation_matrix(theta):
    """2d right-multiply rotation matrix clockwise around origin."""
    ct = np.cos(theta)
    st = np.sin(theta)
    return np.array([[ct, st], [-st, ct]])

# Scale, rotate
samp_points = scale * np.dot(samp_points, right_2d_rotation_matrix(theta)) 
# Translate
samp_points += np.tile(offset, (num_to_keep,1))
# Flip
assert s == 1 or s == -1, s
flip_mat = np.array([[s, 0], [0, 1]])
samp_points = np.dot(samp_points, flip_mat)

print
print 'theta:', theta
print 'scale:', scale
print 'offset:', offset
print 'x parity:', s
print 'samp_points:'
print samp_points
print samp_points.shape

#plt.plot(ref_points[:,0], ref_points[:,1], '.')
#plt.plot(samp_points[:,0], samp_points[:,1], 'r.')
#plt.show()

#plt.plot(ref_samp_points[:3,0], ref_samp_points[:3,1])
#plt.plot(samp_points[:3,0], samp_points[:3,1])
#plt.axis('equal')
#plt.show()

#--------------------------------------------------------------------------------
# Fingerprint ref and samp
#--------------------------------------------------------------------------------

print 'samp_tree...'
samp_tree = scipy.spatial.KDTree(samp_points)

def find_fingerprint(point, source):
    """find_fingerprint(point, source)

    Finds the fingerprint of the given point.
    
    point must be a point from source, where source is either 'ref' or 'samp'.
    
    If a and b are vectors from the point to the first and second closest points, the fingerprint
    is the tuple:
                    (|a|/|b|, np.dot(a,b)/(|a|*|b|))
    That is, the ratio of their lengths and the cosine of the angle between them.
    """
    # Note that the point itself will be the closest, so we use the second and third.
    dists, nbr_idxs = globals().get(source + '_tree').query(point, 3)
    assert dists[0] == 0, dists
    nbrs = globals().get(source + '_points')[nbr_idxs]
    ratio = dists[1]/dists[2]
    angle = np.dot(nbrs[1]-point, nbrs[2]-point) / (dists[1] * dists[2])
    return ratio, angle

def ref_find_fingerprint(x):
    return find_fingerprint(x, 'ref')

def samp_find_fingerprint(x):
    return find_fingerprint(x, 'samp')

pool = mp.Pool(processes = 14)
print 'ref_fing...'
ref_fing = np.array(pool.map(ref_find_fingerprint, ref_points, chunksize=1000))

print 'samp_fing...'
samp_fing = np.array(pool.map(samp_find_fingerprint, samp_points, chunksize=1000))
pool.close()

#--------------------------------------------------------------------------------
# KD Tree of fingerprints
#--------------------------------------------------------------------------------

print 'ref_fing_tree...'
ref_fing_tree = scipy.spatial.KDTree(ref_fing)

def fing_dist_and_idxs(args):
    samp_idx, fing = args
    dist, ref_idx = ref_fing_tree.query(fing)
    return dist, ref_idx, samp_idx

print 'fing_dist_idxs...'
pool = mp.Pool(processes = 14)
fing_dist_idxs = np.array(pool.map(fing_dist_and_idxs, enumerate(samp_fing), chunksize=1000))
pool.close()


sorted_fing_dist_idxs = fing_dist_idxs[np.argsort(fing_dist_idxs[:,0])]

print sorted_fing_dist_idxs

#--------------------------------------------------------------------------------
# Find putative best mappings
#--------------------------------------------------------------------------------

def approx_eq(a, b, tol=1e-3):
    return bool(abs(a-b)<tol)

def ls_mapping_given_all_ref_samp_mappings(ref_samp_mappings):
    """mapping_given_idxs(ref_samp_mappings)

    Input: set of tuples of (ref_idx, samp_idx) mappings.

    Output: scaling lambda, rotation theta, x_offset, y_offset

    We here solve the matrix least squares equation Ax = b, where
                                                                                          
            [ x0r -y0r 1 0 ]
            [ y0r  x0r 0 1 ]
        A = [ x1r -y1r 1 0 ]
            [ y1r  x1r 0 1 ]
                  . . .
            [ xnr -ynr 1 0 ]
            [ ynr  xnr 0 1 ]

    and

        b = [ x0s y0s x1s y1s . . . xns yns ]^T

    The r and s subscripts indicate ref and samp coords.
    
    The interpretation of x is then given by

        x = [ a b x_offset y_offset ]^T

    where
        a = lambda cos(theta), and
        b = lambda sin(theta)

    This system of equations is then finally solved for lambda and theta.
    """
    A = np.zeros((2*len(ref_samp_mappings), 4))
    b = np.zeros((2*len(ref_samp_mappings),))
    for i, (ref_idx, samp_idx) in enumerate(ref_samp_mappings):
        xir, yir = ref_points[ref_idx]
        A[2*i, :]   = [xir, -yir, 1, 0]
        A[2*i+1, :] = [yir,  xir, 0, 1]

        xis, yis = samp_points[samp_idx]
        b[2*i]   = xis
        b[2*i+1] = yis

    x = np.linalg.lstsq(A, b)[0]
    a, b, x_offset, y_offset = x

    theta = np.arctan2(b, a)
    lbda = a / np.cos(theta)

    return lbda, theta, x_offset, y_offset


def ls_mapping_given_idxs(ref_idx, samp_idx):
    """mapping_given_idxs(ref_idx, samp_idx)

    Input: index of ref point and index of corresponding sample point.

    Output: scaling lambda, rotation theta, x_offset, y_offset

    We here solve the matrix least squares equation Ax = b, where
                                                                                          
            [ x0r -y0r 1 0 ]
            [ y0r  x0r 0 1 ]
        A = [ x1r -y1r 1 0 ]
            [ y1r  x1r 0 1 ]
            [ x2r -y2r 1 0 ]
            [ y2r  x2r 0 1 ]

    and

        b = [ x0s y0s x1s y1s x2s y2s ]^T

    The r and s subscripts indicate ref and samp coords.
    
    The interpretation of x is then given by

        x = [ a b x_offset y_offset ]^T

    where
        a = lambda cos(theta), and
        b = lambda sin(theta)

    This system of equations is then finally solved for lambda and theta.
    """
    ref_idx, samp_idx = map(int, (ref_idx, samp_idx))
    if (ref_fing[ref_idx][0] > 0.99
            or abs(ref_fing[ref_idx][1]) > 0.97
            or abs(samp_fing[samp_idx][1]) > 0.97):
        return -1
    ref_dists, ref_nbr_idxs = ref_tree.query(ref_points[ref_idx], 3)
    samp_dists, samp_nbr_idxs = samp_tree.query(samp_points[samp_idx], 3)
    x0r, y0r = ref_points[ref_nbr_idxs[0]]
    x1r, y1r = ref_points[ref_nbr_idxs[1]]
    x2r, y2r = ref_points[ref_nbr_idxs[2]]

    x0s, y0s = samp_points[samp_nbr_idxs[0]]
    x1s, y1s = samp_points[samp_nbr_idxs[1]]
    x2s, y2s = samp_points[samp_nbr_idxs[2]]

    A = np.array([  [ x0r, -y0r, 1, 0 ],
                    [ y0r,  x0r, 0, 1 ],
                    [ x1r, -y1r, 1, 0 ],
                    [ y1r,  x1r, 0, 1 ],
                    [ x2r, -y2r, 1, 0 ],
                    [ y2r,  x2r, 0, 1 ]])
    b = np.array([x0s, y0s, x1s, y1s, x2s, y2s])

    x = np.linalg.lstsq(A, b)[0]
    a, b, x_offset, y_offset = x

    theta = np.arctan2(b, a)
    lbda = a / np.cos(theta)

    return lbda, theta, x_offset, y_offset

def mapping_given_idxs(ref_idx, samp_idx):
    """mapping_given_idxs(ref_idx, samp_idx)

    Input: index of ref point and index of corresponding sample point.

    Output: scaling lambda, rotation theta, x_offset, y_offset

    We here solve the matrix equation Ax = b, where
                                                                                          
            [ x0r  0  -y0r  0  1 0 ]
            [  0  y0r   0  x0r 0 1 ]
        A = [ x1r  0  -y1r  0  1 0 ]
            [  0  y1r   0  x1r 0 1 ]
            [ x2r  0  -y2r  0  1 0 ]
            [  0  y2r   0  x2r 0 1 ]

    and

        b = [ x0s y0s x1s y1s x2s y2s ]^T

    The r and s subscripts indicate ref and samp coords.
    
    The interpretation of x is then given by

        x = [ a aa b bb x_offset y_offset ]^T

    where
        a = aa = lambda cos(theta), and
        b = bb = lambda sin(theta)

    This system of equations is then finally solved for lambda and theta.
    """
    ref_idx, samp_idx = map(int, (ref_idx, samp_idx))
    if (ref_fing[ref_idx][0] > 0.99
            or abs(ref_fing[ref_idx][1]) > 0.97
            or abs(samp_fing[samp_idx][1]) > 0.97):
        return -1
    ref_dists, ref_nbr_idxs = ref_tree.query(ref_points[ref_idx], 3)
    samp_dists, samp_nbr_idxs = samp_tree.query(samp_points[samp_idx], 3)
    x0r, y0r = ref_points[ref_nbr_idxs[0]]
    x1r, y1r = ref_points[ref_nbr_idxs[1]]
    x2r, y2r = ref_points[ref_nbr_idxs[2]]

    x0s, y0s = samp_points[samp_nbr_idxs[0]]
    x1s, y1s = samp_points[samp_nbr_idxs[1]]
    x2s, y2s = samp_points[samp_nbr_idxs[2]]

    A = np.array([  [ x0r,  0,  -y0r,  0,  1, 0 ],
                    [  0,  y0r,   0,  x0r, 0, 1 ],
                    [ x1r,  0,  -y1r,  0,  1, 0 ],
                    [  0,  y1r,   0,  x1r, 0, 1 ],
                    [ x2r,  0,  -y2r,  0,  1, 0 ],
                    [  0,  y2r,   0,  x2r, 0, 1 ]])
    b = np.array([x0s, y0s, x1s, y1s, x2s, y2s])

    x = np.linalg.solve(A, b)
    a, aa, b, bb, x_offset, y_offset = x

    tol = 0.01 * max(map(abs, [a, aa, b, bb]))

    if not (approx_eq(a, aa, tol=tol) and approx_eq(b, bb, tol=tol)):
        return -1

    a = np.mean([a, aa])
    b = np.mean([b, bb])

    theta = np.arctan2(b, a)
    lbda = a / np.cos(theta)

    return lbda, theta, x_offset, y_offset

num_best_fings = 40000
mappings = np.zeros((num_best_fings, 4))
mappings_ref_samp_idxs = np.zeros((num_best_fings, 2), dtype=np.int)
j = 0
for i, (dist, ref_idx, samp_idx) in enumerate(fing_dist_idxs):
    ref_idx, samp_idx = map(int, (ref_idx, samp_idx))
    x = ls_mapping_given_idxs(ref_idx, samp_idx)
    if x == -1:
        continue
    mappings[j] = x
    mappings_ref_samp_idxs[j] = [ref_idx, samp_idx]
    j += 1
    if j >= num_best_fings:
        break
print '%d tries for %d hits' % (i, num_best_fings)

#--------------------------------------------------------------------------------
# Find cluster based on nearest neighbors
#--------------------------------------------------------------------------------

print 'cluster mappings...'
mapping_tree = scipy.spatial.KDTree(mappings)

def triangle_ref_samp_idxs_given_mapping_idx(map_idx):
    """Return list of (ref_idx, samp_idx) tuples for all three points in triangle."""
    ref_idx, samp_idx = mappings_ref_samp_idxs[map_idx]
    ref_dists, ref_nbr_idxs = ref_tree.query(ref_points[ref_idx], 3)
    samp_dists, samp_nbr_idxs = samp_tree.query(samp_points[samp_idx], 3)
    return [(r_idx, s_idx) for r_idx, s_idx in zip(ref_nbr_idxs, samp_nbr_idxs)]


num_top_maps = 10
for k in range(1,12):
    print '\nk =', k
    def mapping_kth_neighbor_dist(args):
        mapping_idx, mapping = args
        dists, nbrs = mapping_tree.query(mapping, k+1)
        return dists[-1], mapping_idx
    
    pool = mp.Pool(processes = 14)
    kth_nbr_dists = np.array(pool.map(mapping_kth_neighbor_dist, enumerate(mappings), chunksize=1000))
    pool.close()
    
    sorted_kth_nbr_dists = kth_nbr_dists[np.argsort(kth_nbr_dists[:,0])]

    ref_samp_idx_mappings = set()
    for i in range(num_top_maps):
        map_idx = int(sorted_kth_nbr_dists[i,1])
        if i < 5:
            print mappings[map_idx]
        ref_samp_idx_mappings.update(triangle_ref_samp_idxs_given_mapping_idx(map_idx))
    
    print 'mean:', np.mean(mappings[map(int, sorted_kth_nbr_dists[:num_top_maps,1])], axis=0)
    print 'median:', np.median(mappings[map(int, sorted_kth_nbr_dists[:num_top_maps,1])], axis=0)
    print 'lsqr:', ls_mapping_given_all_ref_samp_mappings(ref_samp_idx_mappings)
    plt.plot(-np.log10(sorted_kth_nbr_dists[:,0]))
    plt.title('%f,  %f,  %f,  %f' %
            (ls_mapping_given_all_ref_samp_mappings(ref_samp_idx_mappings)))
    plt.show()
sys.exit()


#--------------------------------------------------------------------------------
# Find best mappings from kernel density
#--------------------------------------------------------------------------------

print 'kde...'
kde = neighbors.KernelDensity(bandwidth=30.0)
kde.fit(mappings)
logp = kde.score(mappings)
best_point = mappings[np.argmax(logp)]
print best_point

#kde = scipy.stats.gaussian_kde(mappings)
#p = kde.evaluate(mappings)
#best_point = mappings[np.argmax(p)]
#print best_point

#--------------------------------------------------------------------------------
# Plot results
#--------------------------------------------------------------------------------

#from matplotlib.patches import Ellipse
#esr1 = Ellipse(xy=g.means_[0][:2], width=g.covars_[0][0], height=g.covars_[0][1], facecolor='none')
alpha = 0.01

def plot_1d_kde(ax, data):
    kde = scipy.stats.gaussian_kde(data)
    xx = np.linspace(min(data), max(data), 200)
    yy = kde.evaluate(xx)
    ax.plot(xx, yy)

# 1d kdes
titles = ['lambda', 'theta', 'x offset', 'y offset']
fig = plt.figure()
for i in range(4):
    ax = fig.add_subplot(2,2,i+1)
    plot_1d_kde(ax, mappings[:,i])
    ax.set_title(titles[i])


fig = plt.figure()
ax = fig.add_subplot(121)
Visualize.enhanced_scatter(mappings[:,0], mappings[:,1], ax, hists_height=.2)
ax.plot([best_point[0]], [best_point[1]], 'r*')

ax = fig.add_subplot(122)
Visualize.enhanced_scatter(mappings[:,2], mappings[:,3], ax, hists_height=.2)
ax.plot([best_point[2]], [best_point[3]], 'r*')

# 2d scatter plots
fig = plt.figure()
ax = fig.add_subplot(221)
ax.scatter(mappings[:,0], mappings[:,1], 'o', alpha=alpha)
Visualize.enhanced_scatter(mappings[:,0], mappings[:,1], ax, hist_height=.2)
ax.plot([best_point[0]], [best_point[1]], 'r*')
#ax.add_artist(esr1)

ax = fig.add_subplot(222)
ax.scatter(mappings[:,2], mappings[:,3], 'o', alpha=alpha)
ax.plot([best_point[2]], [best_point[3]], 'r*')

ax = fig.add_subplot(223)
ax.scatter(mappings[:,0], mappings[:,1], 'o', alpha=alpha)
ax.plot([best_point[0]], [best_point[1]], 'r*')
ax.set_xlim((.09, .11))
ax.set_ylim((1.55, 1.6))

ax = fig.add_subplot(224)
ax.scatter(mappings[:,2], mappings[:,3], 'o', alpha=alpha)
ax.plot([best_point[2]], [best_point[3]], 'r*')
ax.set_xlim((0, 200))
ax.set_ylim((100, 300))
plt.show()
