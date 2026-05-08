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

os.makedirs(tmpDir, exist_ok = True)
## set fixed parameters
refDir = "/ix1/wchen/xiangyu/Ref_data/SCENIC/"
path_to_blacklist=f'{refDir}hg38-blacklist.v2.bed.gz'
os.environ['MALLET_MEMORY'] = '200G'
mallet_path="/ix1/wchen/xiangyu/Biosoft/Mallet-202108/bin/mallet"

## set file path
count_file=f'{dataDir}peaks.mtx'
meta_file=f'{dataDir}metadata.tsv'
region_file=f'{dataDir}features.tsv'
outDir = f'{dataDir}outs/'
os.makedirs(outDir, exist_ok = True)

## read in data
cell_data =  pd.read_csv(meta_file, sep='\t')
cell_data = cell_data.set_index("Unnamed: 0")
region_data =  pd.read_csv(region_file, sep='\t', header = None)
count_matrix = pd.DataFrame(
    mmread(count_file).toarray(), 
    index = region_data.iloc[:,0].tolist(),
    columns = cell_data.index.tolist()
)

## Create cisTopic object and add cell information
cistopic_obj = create_cistopic_object(
    fragment_matrix = count_matrix, 
    path_to_blacklist=path_to_blacklist
)

cistopic_obj.add_cell_data(cell_data)
## Run models
from pycisTopic.lda_models import run_cgs_models_mallet
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics = [10, 20, 30, 40, 50],
    n_cpu = n_thread,
    n_iter = 500,
    random_state = 555,
    alpha = 50,
    alpha_by_topic = True,
    eta = 0.1,
    eta_by_topic = False,
    tmp_path = tmpDir,
    save_path = None,
    mallet_path = mallet_path
)

## Save models
with open(outDir+'Mallet_models_500.pkl', 'wb') as f:
  pickle.dump(models, f)

##
from pycisTopic.lda_models import evaluate_models
## model selection
model = evaluate_models(
    models,
    select_model = None, 
    return_model = True, 
    metrics = ['Arun_2010', 'Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
    plot_metrics = False,
    save= outDir + 'model_selection.pdf'
)

## Add model to cisTopicObject
cistopic_obj.add_LDA_model(model)
## Save cisTopic object
with open(outDir + 'cisTopicObject.pkl', 'wb') as f:
  pickle.dump(cistopic_obj, f)


