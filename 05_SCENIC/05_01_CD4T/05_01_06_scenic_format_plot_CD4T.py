## conda activate /ix1/wchen/xiangyu/conda_env/scenicplus
#############
import os
import mudata
from scenicplus.plotting.dotplot import (heatmap_dotplot, generate_dotplot_df)
from IPython import get_ipython
from scenicplus.scenicplus_class import mudata_to_scenicplus
import pandas as pd
##
os.chdir("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/")
## load data
out_path = '03_output/07_SCENIC/CD4T/scplus_pipeline/Snakemake/'
data_path = "03_output/03_clustering/CD4T/WNN_ADT_RNA/"
scplus_mdata_raw = mudata.read(f'{out_path}/scplusmdata.h5mu')
# load pre-filtered eRegulons
e_regulon_filter = pd.read_csv(f"{out_path}/e_regulon_ff_pos.txt") 
scplus_mdata_raw.uns["e_regulon_filter"] = e_regulon_filter
e_regulon_name_trans = pd.read_csv(f"{out_path}/e_regulon_name_trans.txt") 

# add annotation
scplus_mdata = scplus_mdata_raw.copy()
meta_df = pd.read_csv(f"{data_path}/sc_meta.txt", index_col = 0) 
meta_df.index = meta_df.index + "___cisTopic"
scplus_mdata.mod["scRNA_counts"].obs["ann_level3_final"] = meta_df[["ann_level3_final"]]
scplus_mdata.mod["scRNA_counts"].obs["ann_level4_final"] = meta_df[["ann_level4_final"]]

##
scplus_mdata = mudata.MuData(scplus_mdata.mod)
scplus_mdata.uns = scplus_mdata_raw.uns

## calculate rss
from scenicplus.RSS import (regulon_specificity_scores, plot_rss)
rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scRNA_counts:ann_level4_final",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC",
    "direct_region_based_AUC", "extended_region_based_AUC"],
    selected_regulons = e_regulon_name_trans["signature"].to_list()
)
rss.to_csv(f"{out_path}/rss_CD4_ann_level4_final.txt", header=True, index=True, sep='\t', mode='w')
rssx = plot_rss(
    data_matrix = rss,
    top_n = 5,
    num_columns = 3,
    save = f"{out_path}/rss_ann_level4_final.pdf"
)
## heatmap for gene and atac activity
heatx = heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scRNA_counts:ann_level4_final",
    eRegulon_metadata_key = "e_regulon_filter",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "vertical",
    figsize = (20, 10)
)
heatx.save(f"{out_path}/heat_CD4_ann_level4_final.pdf", format="pdf")

##### Treg subset #####
## subset
scplus_mdata_treg = scplus_mdata[scplus_mdata.mod["scRNA_counts"].obs["ann_level3_final"][:] == "CD4_Treg",:]
scplus_mdata_treg.uns = scplus_mdata.uns
## calculate rss
rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata_treg,
    variable = "scRNA_counts:ann_level4_final",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC",
    "direct_region_based_AUC", "extended_region_based_AUC"],
    selected_regulons = e_regulon_name_trans["signature"].to_list()
)
rss.to_csv(f"{out_path}/rss_Treg_ann_level4_final.txt", header=True, index=True, sep='\t', mode='w')
rssx = plot_rss(
    data_matrix = rss,
    top_n = 5,
    num_columns = 3,
    save = f"{out_path}/rss_Treg_ann_level4_final.pdf"
)

## heatmap for gene and atac activity
size_matrix = scplus_mdata_treg["direct_region_based_AUC"].to_df()
color_matrix = scplus_mdata_treg["direct_gene_based_AUC"].to_df()
group_by = scplus_mdata_treg.obs["scRNA_counts:ann_level4_final"].tolist()
size_features, color_features, feature_names = scplus_mdata_treg.uns["e_regulon_filter"][["Region_signature_name", "Gene_signature_name", "eRegulon_name"]] \
    .drop_duplicates().values.T
#
plotting_df = generate_dotplot_df(
        size_matrix = size_matrix,
        color_matrix = color_matrix,
        group_by = group_by,
        size_features = size_features,
        color_features = color_features,
        feature_names = feature_names,
        scale_size_matrix = True,    
        scale_color_matrix = True,
        group_name = "scRNA_counts:ann_level4_final",
        size_name = "direct_region_based_AUC",
        color_name = "direct_gene_based_AUC",
        feature_name = "eRegulon_name"
)
plotting_df.to_csv(f"{out_path}/heatdot_Treg_ann_level4_final.txt", header=True, index=True, sep='\t', mode='w')
heatx = heatmap_dotplot(
    scplus_mudata = scplus_mdata_treg,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scRNA_counts:ann_level4_final",
    eRegulon_metadata_key = "e_regulon_filter",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "vertical",
    figsize = (6, 20)
)
heatx.save(f"{out_path}/heat_Treg_ann_level4_final.pdf", format="pdf")

##### CD103 TRM subset #####
scplus_mdata_cd103trm = scplus_mdata[scplus_mdata.mod["scRNA_counts"].obs["ann_level4_final"][:] == "CD4_CD103_TRM_Th17",:]
data_path2 = "03_output/03_clustering/CD4_CD103_TRM/WNN_RNA_ATAC/"
meta_df_cd103trm = pd.read_csv(f"{data_path2}/sc_meta.txt", index_col = 0) 
meta_df_cd103trm.index = meta_df_cd103trm.index + "___cisTopic"
meta_df_cd103trm = meta_df_cd103trm.loc[scplus_mdata_cd103trm.obs.index,:]
scplus_mdata_cd103trm.mod["scRNA_counts"].obs["ann_level5_final"] = meta_df_cd103trm[["ann_level5_final"]]
scplus_mdata_cd103trm = mudata.MuData(scplus_mdata_cd103trm.mod)
scplus_mdata_cd103trm.uns = scplus_mdata.uns
## calculate rss
rss_cd103trm = regulon_specificity_scores(
    scplus_mudata = scplus_mdata_cd103trm,
    variable = "scRNA_counts:ann_level5_final",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC",
    "direct_region_based_AUC", "extended_region_based_AUC"]
)
rss_cd103trm.to_csv(f"{out_path}/rss_cd103trm17_ann_level5_final.txt", header=True, index=True, sep='\t', mode='w')
rssx_cd103trm = plot_rss(
    data_matrix = rss_cd103trm,
    top_n = 5,
    num_columns = 2,
    save = f"{out_path}/rss_cd103trm17_ann_level5_final.pdf"
)
## heatmap for gene and atac activity
size_matrix = scplus_mdata_cd103trm["direct_region_based_AUC"].to_df()
color_matrix = scplus_mdata_cd103trm["direct_gene_based_AUC"].to_df()
group_by = scplus_mdata_cd103trm.obs["scRNA_counts:ann_level5_final"].tolist()
size_features, color_features, feature_names = scplus_mdata_cd103trm.uns["e_regulon_filter"][["Region_signature_name", "Gene_signature_name", "eRegulon_name"]] \
    .drop_duplicates().values.T
#
plotting_df = generate_dotplot_df(
        size_matrix = size_matrix,
        color_matrix = color_matrix,
        group_by = group_by,
        size_features = size_features,
        color_features = color_features,
        feature_names = feature_names,
        scale_size_matrix = True,    
        scale_color_matrix = True,
        group_name = "scRNA_counts:ann_level5_final",
        size_name = "direct_region_based_AUC",
        color_name = "direct_gene_based_AUC",
        feature_name = "eRegulon_name"
)
plotting_df.to_csv(f"{out_path}/heatdot_CD103trm17_ann_level5_final.txt", header=True, index=True, sep='\t', mode='w')
#
heatx = heatmap_dotplot(
    scplus_mudata = scplus_mdata_cd103trm,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scRNA_counts:ann_level5_final",
    eRegulon_metadata_key = "e_regulon_filter",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "vertical",
    figsize = (15, 20)
)
heatx.save(f"{out_path}/heat_cd103trm17_ann_level5_final.pdf", format="pdf")

