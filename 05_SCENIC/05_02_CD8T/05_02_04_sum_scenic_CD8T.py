## conda activate /ix1/wchen/xiangyu/conda_env/scenicplus
#############
import os
import mudata
from scenicplus.plotting.dotplot import heatmap_dotplot
from IPython import get_ipython
from scenicplus.scenicplus_class import mudata_to_scenicplus
import pandas as pd

os.chdir("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/")
ct = "CD8T"
out_path = f'03_output/07_SCENIC/{ct}/scplus_pipeline/Snakemake/'
scplus_mdata = mudata.read(f'{out_path}/scplusmdata.h5mu')

# Simplifying and filtering SCENIC+ output
from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
apply_std_filtering_to_eRegulons(scplus_mdata)
e_regulon_metadata_filtered = scplus_mdata.uns["e_regulon_metadata_filtered"]
e_regulon_metadata_filtered.to_csv(f"{out_path}/e_regulon_metadata_filtered.txt", 
                                   header=True, index=False, sep='\t', mode='w')

