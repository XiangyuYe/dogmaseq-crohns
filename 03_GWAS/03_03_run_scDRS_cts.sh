## 0. set parameters
grch=37.3
eth=eur
##
project_path=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/
out_path=${project_path}03_output/01_scDRS/
##
magma=/ix1/wchen/xiangyu/Biosoft/MAGMA/magma
ref_path=/ix1/wchen/xiangyu/Ref_data/
sc_path=${project_path}/03_output/01_scDRS/00_scRNA/

##
p=IBD
n_top=1000
n_ctrl=1000

gs_file=${out_path}/02_gs/${p}_top${n_top}.gs
sc_file=${sc_path}/all/count.h5ad
integ_path=${out_path}/03_integrate/
dstream_path=${out_path}/04_downstream/

##############################################################
module load r/4.5.0
#
integ_pathx=${integ_path}/
dstream_pathx=${dstream_path}/
cov_filex=${sc_path}/all/cov.tsv
mkdir -p ${integ_pathx}
mkdir -p ${dstream_pathx}

# compute-score
scdrs compute-score \
--h5ad-file ${sc_file} \
--h5ad-species human \
--gs-file ${gs_file} \
--gs-species human \
--out-folder ${integ_pathx} \
--cov-file ${cov_filex} \
--flag-filter-data True \
--flag-raw-count True \
--n-ctrl ${n_ctrl} \
--flag-return-ctrl-raw-score False \
--flag-return-ctrl-norm-score True

# sample_ct_level2 sample_ct_level3 sample_ct_level4 ann_level2_refine ann_level3_final ann_level4_final
for ann_col in sample_ct_level2 sample_ct_level4 ann_level2_refine ann_level4_final
do
# compare across celltypes
scdrs perform-downstream \
--h5ad-file ${sc_file} \
--score-file ${integ_pathx}/${p}.full_score.gz \
--out-folder ${dstream_pathx} \
--group-analysis ${ann_col} \
--flag-filter-data True \
--flag-raw-count True

done

# module unload r/4.5.0


