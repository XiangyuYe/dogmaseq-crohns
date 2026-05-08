## conda activate microbio
import scanpy as sc
import pandas as pd
import numpy as np
from cnmf import cNMF, Preprocess
#Initialize the Preprocess object
p = Preprocess(random_seed = 20250528)
# load data
proj_path = "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
ct = "CD4_CD103_TRM"
data_path = f"{proj_path}03_output/03_clustering/{ct}/cNMF/"
##
adata = sc.read_10x_mtx(
    f'{data_path}/gex/',
    var_names = 'gene_symbols',
    cache = False)
meta_df = pd.read_csv(f"{data_path}/metadata.tsv", index_col = 0) 
adata.obs["Batch"] = meta_df[["Batch"]]
np.random.seed(20250528)

# Batch correct the data and save the corrected high-variance gene data to adata_c, and the TPM normalized data to adata_tpm 
(adata_c, adata_tpm, hvgs) = p.preprocess_for_cnmf(
                             adata, 
                             harmony_vars = ["Batch"], 
                             n_top_rna_genes = 3000, 
                             save_output_base = f"{data_path}/{ct}_harmony")

#Then run cNMF passing in the corrected counts file, tpm_fn, and HVGs as inputs
cnmf_obj_corrected = cNMF(output_dir = f"{data_path}/", name = f'BatchCorrected_{ct}_top3000')
cnmf_obj_corrected.prepare(counts_fn = f'{data_path}/{ct}_harmony.Corrected.HVG.Varnorm.h5ad',
                           tpm_fn = f'{data_path}/{ct}_harmony.TP10K.h5ad',
                           genes_file = f'{data_path}/{ct}_harmony.Corrected.HVGs.txt',
                           components = np.arange(5,20), 
                           n_iter = 100, 
                           seed = 20250528, 
                           num_highvar_genes = 3000)

#Then proceed with the rest of cNMF as normal
cnmf_obj_corrected.factorize(worker_i = 0, total_workers = 1)
cnmf_obj_corrected.combine()
cnmf_obj_corrected.k_selection_plot(close_fig = True)

##### determine for CD4_CD103_TRM #####
# k1=9
k_sel = 11
cnmf_obj_corrected.consensus(k = k_sel, density_threshold = 2, build_ref = True)
dens_thresh = 0.02
cnmf_obj_corrected.consensus(k = k_sel, density_threshold = dens_thresh, build_ref = True)
usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj_corrected.load_results(K = k_sel, density_threshold = dens_thresh)

# output
col_char = ["GEP" + str(num) for num in range(1, k_sel + 1)]
usage.columns = col_char
spectra_scores.columns = col_char
spectra_tpm.columns = col_char
top_genes.columns = col_char
usage.to_csv(f'{data_path}/{ct}_harmony_usage_k{k_sel}.txt', header = True, index = True, sep = '\t', mode = 'w')
spectra_scores.to_csv(f'{data_path}/{ct}_harmony_spectra_scores_k{k_sel}.txt', header = True, index = True, sep = '\t', mode = 'w')
spectra_tpm.to_csv(f'{data_path}/{ct}_harmony_spectra_tpm_k{k_sel}.txt', header = True, index = True, sep = '\t', mode = 'w')
top_genes.to_csv(f'{data_path}/{ct}_harmony_top_genes_k{k_sel}.txt', header = True, index = False, sep = '\t', mode = 'w')

##### determine for CD4_CD103_TRM #####
k_sel = 16
cnmf_obj_corrected.consensus(k = k_sel, density_threshold = 2, build_ref = True)
dens_thresh = 0.1
cnmf_obj_corrected.consensus(k = k_sel, density_threshold = dens_thresh, build_ref = True)
usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj_corrected.load_results(K = k_sel, density_threshold = dens_thresh)

# output
col_char = ["GEP" + str(num) for num in range(1, k_sel + 1)]
usage.columns = col_char
spectra_scores.columns = col_char
spectra_tpm.columns = col_char
top_genes.columns = col_char
usage.to_csv(f'{data_path}/{ct}_harmony_usage_k{k_sel}.txt', header = True, index = True, sep = '\t', mode = 'w')
spectra_scores.to_csv(f'{data_path}/{ct}_harmony_spectra_scores_k{k_sel}.txt', header = True, index = True, sep = '\t', mode = 'w')
spectra_tpm.to_csv(f'{data_path}/{ct}_harmony_spectra_tpm_k{k_sel}.txt', header = True, index = True, sep = '\t', mode = 'w')
top_genes.to_csv(f'{data_path}/{ct}_harmony_top_genes_k{k_sel}.txt', header = True, index = False, sep = '\t', mode = 'w')

##
k_sel = 8
cnmf_obj_corrected.consensus(k = k_sel, density_threshold = 2, build_ref = True)
dens_thresh = 0.1
cnmf_obj_corrected.consensus(k = k_sel, density_threshold = dens_thresh, build_ref = True)
usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj_corrected.load_results(K = k_sel, density_threshold = dens_thresh)

# output
col_char = ["GEP" + str(num) for num in range(1, k_sel + 1)]
usage.columns = col_char
spectra_scores.columns = col_char
spectra_tpm.columns = col_char
top_genes.columns = col_char
usage.to_csv(f'{data_path}/{ct}_harmony_usage_k{k_sel}.txt', header = True, index = True, sep = '\t', mode = 'w')
spectra_scores.to_csv(f'{data_path}/{ct}_harmony_spectra_scores_k{k_sel}.txt', header = True, index = True, sep = '\t', mode = 'w')
spectra_tpm.to_csv(f'{data_path}/{ct}_harmony_spectra_tpm_k{k_sel}.txt', header = True, index = True, sep = '\t', mode = 'w')
top_genes.to_csv(f'{data_path}/{ct}_harmony_top_genes_k{k_sel}.txt', header = True, index = False, sep = '\t', mode = 'w')
