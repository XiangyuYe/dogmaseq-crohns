import scanpy as sc
import pandas as pd

# Load default TCAT reference from starCAT databse
# tcat = starCAT(reference='TCAT.V1')
# Load cell x genes counts data
data_path = "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/03_output/01_scDRS/00_scRNA/"
####### all #######
adata = sc.read_10x_mtx(
    f'{data_path}all/',
    var_names = 'gene_symbols',
    cache = False)
#
df = pd.read_csv(f'{data_path}all/metadata.tsv', index_col=0)
df = df.loc[adata.obs_names]
adata.obs = df
##
adata.write(f'{data_path}all/count.h5ad')
####### CD4T #######
adata = sc.read_10x_mtx(
    f'{data_path}CD4/',
    var_names = 'gene_symbols',
    cache = False)
#
df = pd.read_csv(f'{data_path}CD4/metadata.tsv', index_col=0)
df = df.loc[adata.obs_names]
adata.obs = df
##
adata.write(f'{data_path}CD4/count.h5ad')
####### CD4_CD103_TRM #######
adata = sc.read_10x_mtx(
    f'{data_path}CD4_CD103_TRM/',
    var_names = 'gene_symbols',
    cache = False)
#
df = pd.read_csv(f'{data_path}CD4_CD103_TRM/metadata.tsv', index_col=0)
df = df.loc[adata.obs_names]
adata.obs = df
##
adata.write(f'{data_path}CD4_CD103_TRM/count.h5ad')
####### CD4_TRM #######
adata = sc.read_10x_mtx(
    f'{data_path}CD4_TRM/',
    var_names = 'gene_symbols',
    cache = False)
#
df = pd.read_csv(f'{data_path}CD4_TRM/metadata.tsv', index_col=0)
df = df.loc[adata.obs_names]
adata.obs = df
##
adata.write(f'{data_path}CD4_TRM/count.h5ad')
