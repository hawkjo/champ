"""
Variables and methods for finding local files and folders.
"""

import os
import numpy as np

base_dir = '/home/hawkjo/BillProjects/constellation'
fig_dir = os.path.join(base_dir, 'figs')
jah_base_dir = '/home/jah/projects/ilya/experiments'
jah_04_26_dir = os.path.join(jah_base_dir, '2015_04_26_imaging_data')

dname_given_project_name = {
        'SA15064': '150331_NS500358_0049_AH75YLBGXX',
        'SA15060': '150401_M02288_0122_000000000-AC6PA',
        'SA15063': '150401_NS500358_0050_AH77JCBGXX',
        'SA15067': '150403_M01012_0037_000000000-AELH5',
        'SA15066': '150407_M01012_0038_000000000-AEU9J',
        'SA15069': '150407_NS500358_0051_AH75Y3BGXX',
        'SA15065': '150408_NS500358_0052_AH77TGBGXX',
        'SA15070': '150410_NS500358_0053_AH773VBGXX',
        }

def phiX_rcs_given_project_name(pname):
    fpath = os.path.join(jah_base_dir,
                         dname_given_project_name[pname],
                         'phiX_mappings/xys.npz')
    return np.load(fpath)
