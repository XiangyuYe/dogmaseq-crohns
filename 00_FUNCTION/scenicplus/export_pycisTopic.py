## never name your script the same as some module!!! (pycisTopic.py to run_pycisTopic.py)
import warnings
warnings.simplefilter(action='ignore')
import pycisTopic
import os
from pycisTopic.cistopic_class import *
from scipy.io import mmread
import pickle 
import sys
import getopt
import pandas as pd

long_opts_list = ['dataDir=', 'tmpDir=', 'anno_col=', 'n_thread=']
param_dict = {'dataDir': None, 'tmpDir': None, 'anno_col': None, 'n_thread': None}
opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)

for opt, arg in opts:
    if opt == "--dataDir": 
        param_dict['dataDir'] = arg
    elif opt == "--tmpDir": 
        param_dict['tmpDir'] = arg
    elif opt == "--anno_col": 
        param_dict['anno_col'] = arg  
    elif opt == "--n_thread": 
        param_dict['n_thread'] = arg          
# set parameters
dataDir = param_dict['dataDir']
tmpDir = param_dict['tmpDir']
anno_col = param_dict['anno_col']
n_thread = int(param_dict['n_thread'])

## input parameters 
# dataDir = "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/03_output/07_SCENIC/CD4/small_intestine/"
# tmpDir = "/ix1/wchen/xiangyu/Tmp/PC/small_intestine/"
# tmpDir2 = "/ix1/wchen/xiangyu/Tmp/"
# anno_col = "ann2"
# n_thread = 20

outDir = f'{dataDir}outs/'
os.makedirs(outDir, exist_ok = True)
os.makedirs(tmpDir, exist_ok = True)
## Topic binarization & QC
# Load cisTopic object
import pickle
infile = open(outDir + 'cisTopicObject.pkl', 'rb')
cistopic_obj = pickle.load(infile)
infile.close()
## Topic binarization & QC
from pycisTopic.topic_binarization import binarize_topics
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, 
    method='ntop', 
    ntop = 3000,
    plot=False, 
    num_columns=5
)

#
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, 
    method='otsu',
    plot=False, 
    num_columns=5
)

#
binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=False,
    num_columns=5, 
    nbins=100
)

## Differentially Accessible Regions (DARs)
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features)
import numpy as np
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)

#
normalized_imputed_acc_obj = normalize_scores(
    imputed_acc_obj, 
    scale_factor=10**4
)

#
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=False
)

# #
# markers_dict= find_diff_features(
    # cistopic_obj,
    # imputed_acc_obj,
    # variable=anno_col,
    # var_features=variable_regions,
    # contrasts=None,
    # adjpval_thr=0.05,
    # log2fc_thr=np.log2(1.5),
    # n_cpu=n_thread,
    # _temp_dir=tmpDir,
    # split_pattern = '-'
# )

## Save region sets
os.makedirs(os.path.join(outDir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(outDir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(outDir, "region_sets", "Topics_top_3k"), exist_ok = True)
# os.makedirs(os.path.join(outDir, "region_sets", "DARs_cell_type"), exist_ok = True)
#
from pycisTopic.utils import region_names_to_coordinates
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(outDir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
)

#
for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(outDir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
)

# #
# for cell_type in markers_dict:
    # region_names_to_coordinates(
        # markers_dict[cell_type].index
    # ).sort_values(
        # ["Chromosome", "Start", "End"]
    # ).to_csv(
        # os.path.join(outDir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        # sep = "\t",
        # header = False, index = False
# )



