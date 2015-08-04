import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import sklearn.lda
import sklearn.qda
import sklearn.neighbors
import sklearn.tree
import sklearn.svm
import sklearn.ensemble
import sklearn.cross_validation
import sklearn.metrics
import tifffile
import random
import local_config
import image_processing_tools
import multiprocessing as mp
from pprint import pprint
from pathos.multiprocessing import ProcessingPool


#------------------------------------------------------------------------------------------
# Loading data
#------------------------------------------------------------------------------------------

def result_and_image_fpaths(results_dir, im_dir, bname_fmt, bad_images=[], verbose=True):
    """Returns list of result/image tuples. Optional list of bad image ints for removal."""
    bad_results_fnames = ['{0}_all_read_rcs.txt'.format(bname_fmt % i) for i in bad_images]
    results_fpaths = [fpath for fpath in glob.glob(os.path.join(results_dir, '*_all_read_rcs.txt'))
                      if os.path.basename(fpath) not in bad_results_fnames]
    if verbose:
        print 'Num results files:', len(results_fpaths)
        
    res_im_fpaths = []
    for res_fpath in results_fpaths:
        fname = os.path.basename(res_fpath)
        bname = fname[:fname.index('_')]
        res_im_fpaths.append((res_fpath, os.path.join(im_dir, bname + '.tif')))
    return res_im_fpaths


def get_phiX_read_names(project_name):
    fastq_tiles = local_config.phiX_read_names_given_project_name(project_name)
    phiX_read_names = set()
    for tile_key, read_names in fastq_tiles.items():
        phiX_read_names.update(read_names)
    return phiX_read_names


#------------------------------------------------------------------------------------------
# Useful training functions
#------------------------------------------------------------------------------------------
clf_funcs = {
    'lda': sklearn.lda.LDA,
    'qda': sklearn.qda.QDA,
    'knn': sklearn.neighbors.KNeighborsClassifier,
    'svm': sklearn.svm.SVC,
    'dt': sklearn.tree.DecisionTreeClassifier,
    'adaboost': sklearn.ensemble.AdaBoostClassifier,
}


def make_train_classifier(classifier_type, classifiers_dict):
    # Given a dict of classifiers and their functions/kwargs, as constructed here, create a poolable function.
    clf_func = clf_funcs[classifier_type]
    def train_classifier(c_key):
        d = classifiers_dict[c_key]
        clf = clf_func(**d['kwargs'])
        clf.fit(Xtrain, ytrain)
        ypred = clf.predict(Xtest)
        yscore = getattr(clf, d['score_func'])(Xtest)
        return clf, ypred, yscore
    return train_classifier

def make_train_classifier_results(trn_clf_func):
    def train_classifier_results(c_key):
        clf, ypred, yscore = trn_clf_func(c_key)
        results = (
            c_key,
            sklearn.metrics.accuracy_score(ytest, ypred),
            sklearn.metrics.precision_score(ytest, ypred),
            sklearn.metrics.recall_score(ytest, ypred),
            sklearn.metrics.f1_score(ytest, ypred),
        )
        if c_key[0] not in ['KNN', 'DT']:
            results = results + (sklearn.metrics.roc_auc_score(ytest, yscore),)
        else:
            results = results + (None,)
        return results
    return train_classifier_results

def run_all_classifiers(key, classifiers):
    train_classifier_results = make_train_classifier_results(make_train_classifier(key, classifiers))
    pool = ProcessingPool(12)
    return pool.map(train_classifier_results, classifiers.keys())


def print_results(results, num_top=5):
    header = '%40s   %5s   %5s   %5s   %5s   %5s' % ('Run', 'Acc', 'Prec', 'Rcl', 'F1', 'AUC')
    def fmt_n(n):
        if n is None:
            return '%5s' % None
        else:
            return '%.3f' % n
    print 'Best by F1-Score:'
    print
    results.sort(key=lambda tup: -tup[4])
    print header
    for run in results[:num_top]:
        print '   '.join(['%40s' % (run[0],)] + [fmt_n(n) for n in run[1:]])
    if not set([None]) == set([run[5] for run in results]):
        print
        print 'Best by AUC:'
        print
        print header
        results.sort(key=lambda tup: tup[5], reverse=True)
        for run in results[:num_top]:
            print '   '.join(['%40s' % (run[0],)] + [fmt_n(n) for n in run[1:]])


def stoftoi(s):
    return int(float(s))
    

def get_train_test_data(
        res_im_fpaths,
        on_read_names,
        side_px=1,
        include_dist_to_center=True,
        include_rc=True,
        filter_image=False,
        verbose=False,
        samp_size=20000,
        test_size=0.3,
        ):
    X = []
    y = []
    for i, (res_fpath, im_fpath) in enumerate(res_im_fpaths):
        if verbose:
            print i,
            sys.stdout.flush()
        im = tifffile.imread(im_fpath)
        med = np.median(im)
        im = im / float(med)
        im -= 1
        if filter_image:
            filt = image_processing_tools.signal_hist_and_func(im, verbose=False, plot_curves=False)
            im = im * filt(im)
            
        def dist_to_im_center(r, c):
            return np.linalg.norm(np.array(im.shape)/2 - np.array([r, c]))
        
        for line in open(res_fpath):
            read_name, r, c = line.strip().split()
            r, c = map(stoftoi, (r, c))
            if side_px <= r < im.shape[0] - side_px - 1 and side_px <= c < im.shape[0] - side_px - 1:
                xx = im[r-side_px:r+side_px+1, c-side_px:c+side_px+1].flatten()
                if include_dist_to_center:
                    xx = np.r_[xx, dist_to_im_center(r, c)]
                if include_rc:
                    xx = np.r_[xx, r, c]
                X.append(xx)
                y.append(int(read_name in phiX_read_names))
                
    X = np.array(X)
    y = np.array(y)
    
    sub_idxs = random.sample(range(len(y)), samp_size)
    subX = X[sub_idxs]
    suby = y[sub_idxs]
    Xtrain, Xtest, ytrain, ytest = sklearn.cross_validation.train_test_split(
            subX, suby, test_size=test_size, random_state=42)
    if verbose:
        print
        print 'Data shapes'
        print 'X and y:', X.shape, y.shape
        print 'Xtrain and ytrain:', Xtrain.shape, ytrain.shape
        print 'Xtest and ytest:', Xtest.shape, ytest.shape
