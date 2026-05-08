# conda activate /ix1/wchen/xiangyu/conda_env/starCAT
import pandas as pd
import numpy as np
import scanpy as sc
import os
import matplotlib.pyplot as plt
import sys
import seaborn as sns
from starcat import starCAT
#
proj_path = "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
# Load default TCAT reference from starCAT databse
tcat = starCAT(reference = 'TCAT.V1')
# Load cell x genes counts data
data_path = f"{proj_path}03_output/02_clean/all/CD4T/RNA/"
out_path = f"{proj_path}03_output/03_clustering/CD4T/TCAT/"
os.makedirs(out_path, exist_ok = True)

adata = sc.read_10x_mtx(
    f"{data_path}",
    var_names = 'gene_symbols',
    cache = False)

# Run starCAT
usage, scores = tcat.fit_transform(adata)
# output
usage.to_csv(f"{out_path}usage.txt", header = True, index = True, sep = '\t', mode = 'w')
scores.to_csv(f"{out_path}scores.txt", header = True, index = True, sep = '\t', mode = 'w')
