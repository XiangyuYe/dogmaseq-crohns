## 0. set parameters
grch=37.3
eth=eur
#
project_path=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/
out_path=${project_path}03_output/01_scDRS/
#
magma=/ix1/wchen/xiangyu/Biosoft/MAGMA/magma
zscore_format=${project_path}code/00_scDRS_format_zscore.R
ref_path=/ix1/wchen/xiangyu/Ref_data/

mkdir -p ${out_path}/01_magma
mkdir -p ${out_path}/02_gs
mkdir -p ${out_path}/03_integrate
mkdir -p ${out_path}/04_downstream

module load gcc/12.2.0
p=CD
n_top=1000
n_ctrl=1000
summ_file=${ref_path}GWAS/GRCh37/${p}.txt
sc_path=${project_path}/03_output/01_scDRS/00_scRNA/
##############################################################

for ann_col in sample_ct_level4 ann_level4_refine
do

# compare across celltypes
scdrs perform-downstream \
--h5ad-file ${sc_path}/CD4/count.h5ad \
--score-file ${out_path}/03_integrate/CD4T/all_sub/${covx}/${p}.full_score.gz \
--out-folder ${out_path}/04_downstream/CD4T/all_sub/${covx}/ \
--group-analysis ${ann_col} \
--gene-analysis \
--flag-filter-data True \
--flag-raw-count True

done


################################################################
covx=cov1

# compute-score
sc_path=${project_path}/03_output/01_scDRS/00_scRNA/
scdrs compute-score \
--h5ad-file ${sc_path}/CD4_CD103_TRM/count.h5ad \
--h5ad-species human \
--gs-file ${out_path}/02_gs/${p}_top${n_top}.gs \
--gs-species human \
--out-folder ${out_path}/03_integrate/CD4_CD103_TRM/${covx}/ \
--cov-file ${sc_path}/CD4_CD103_TRM/${covx}.tsv \
--flag-filter-data True \
--flag-raw-count True \
--n-ctrl ${n_ctrl} \
--flag-return-ctrl-raw-score False \
--flag-return-ctrl-norm-score True

# # compute-score
# sc_path=${project_path}/03_output/01_scDRS/00_scRNA/
# scdrs compute-score \
# --h5ad-file ${sc_path}/CD4_TRM/count.h5ad \
# --h5ad-species human \
# --gs-file ${out_path}/02_gs/${p}_top${n_top}.gs \
# --gs-species human \
# --out-folder ${out_path}/03_integrate/CD4_TRM/${covx}/ \
# --cov-file ${sc_path}/CD4_TRM/${covx}.tsv \
# --flag-filter-data True \
# --flag-raw-count True \
# --n-ctrl ${n_ctrl} \
# --flag-return-ctrl-raw-score False \
# --flag-return-ctrl-norm-score True

for ann_col in ann_level4_refine sample_ct_level4
do
# compare across celltypes
scdrs perform-downstream \
--h5ad-file ${sc_path}/CD4_CD103_TRM/count.h5ad \
--score-file ${out_path}/03_integrate/CD4_CD103_TRM/all_sub/${covx}/${p}.full_score.gz \
--out-folder ${out_path}/04_downstream/CD4_CD103_TRM/all_sub/${covx}/ \
--group-analysis ${ann_col} \
--gene-analysis \
--flag-filter-data True \
--flag-raw-count True


# compare across celltypes
scdrs perform-downstream \
--h5ad-file ${sc_path}/CD4_TRM/count.h5ad \
--score-file ${out_path}/03_integrate/CD4_CD103_TRM/CD4T_sub/${covx}/${p}.full_score.gz \
--out-folder ${out_path}/04_downstream/CD4_CD103_TRM/CD4T_sub/${covx}/ \
--group-analysis ${ann_col} \
--gene-analysis \
--flag-filter-data True \
--flag-raw-count True

done
