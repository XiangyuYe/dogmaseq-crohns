import scanpy as sc
import pandas as pd
##
ct="CD4T"
dataDir = f'/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/03_output/07_SCENIC/{ct}/'
countDir=f'{dataDir}/gex_count/'
meta_file=f'{dataDir}metadata.tsv'
adata = sc.read_10x_mtx(
    countDir,
    var_names = "gene_symbols"
)
adata.var_names_make_unique()
cell_data =  pd.read_csv(meta_file, sep='\t')
cell_data = cell_data.set_index("Unnamed: 0")
adata.obs = cell_data.loc[adata.obs_names]
## normalize data (WITHOUT REMOVING RAW COUNTS) 
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write(f'{dataDir}adata.h5ad')

