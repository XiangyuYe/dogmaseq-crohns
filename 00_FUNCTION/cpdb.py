##
import sys
import os
import getopt

long_opts_list = ['counts_file_path=', 'meta_file_path=', 'active_tf_path=', 'out_file_path=', 'nthread=']
param_dict = {'counts_file_path': None, 'meta_file_path': None, 'active_tf_path': None, 'out_file_path': None, 'nthread': 1}
opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)

for opt, arg in opts:
    if opt == "--counts_file_path": 
        param_dict['counts_file_path'] = arg
    elif opt == "--meta_file_path": 
        param_dict['meta_file_path'] = arg
    elif opt == "--active_tf_path": 
        param_dict['active_tf_path'] = arg
    elif opt == "--out_file_path": 
        param_dict['out_file_path'] = arg
    elif opt == "--nthread": 
        param_dict['nthread'] = arg
# set parameters
counts_file_path = param_dict['counts_file_path']
meta_file_path = param_dict['meta_file_path']
out_file_path = param_dict['out_file_path']
active_tf_path = param_dict['active_tf_path']
nthread = int(param_dict['nthread'])
# counts_file_path = '/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/03_output/05_Interaction/cpdb/Sample/'
# out_path = '/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/03_output/05_Interaction/cpdb/Sample/'
# itemx = 'IBD-009_CD_I_I_TI'
# nthread = 2
print(counts_file_path)
print(meta_file_path)
print(out_file_path)
print(active_tf_path)
print(nthread)
# set input
cpdb_file_path = '/ix1/wchen/xiangyu/Biosoft/CellphoneDB/NatureProtocols2024_case_studies/v5.0.0/cellphonedb.zip'
if not os.path.exists(out_file_path):
    os.mkdir(out_file_path)

## run cpdb
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
    active_tfs_file_path = active_tf_path,           # optional: defines cell types and their active TFs.
    microenvs_file_path = None,                      # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                       # optional: whether to score interactions or not. 
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = nthread,                               # number of threads to use in the analysis.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    output_path = out_file_path,                     # Path to save results.
    output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )


