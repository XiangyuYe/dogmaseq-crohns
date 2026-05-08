# conda activate /ix1/wchen/xiangyu/conda_env/scenicplus
## 0. set path
proj_path=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/
ctx=CD8T
## 1. pycisTopic
run_pycisTopic=${proj_path}code/FUNCTION/scenicplus/run_pycisTopic.py
cd /ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/03_output/07_SCENIC/${ctx}/
python ${run_pycisTopic} \
--dataDir ${proj_path}03_output/07_SCENIC/${ctx}/ \
--tmpDir ${proj_path}03_output/07_SCENIC/${ctx}/tmp/ \
--anno_col "ann_raw" \
--n_thread 20
## 2. export pycisTopic
export_pycisTopic=${proj_path}code/FUNCTION/scenicplus/export_pycisTopic.py
python ${export_pycisTopic} \
--dataDir ${proj_path}03_output/07_SCENIC/${ctx}/ \
--tmpDir ${proj_path}03_output/07_SCENIC/${ctx}/tmp/ \
--anno_col "ann_raw" \
--n_thread 20

#### 3. snakemake ####
ct=CD4T
dataDir=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/03_output/07_SCENIC/${ct}/

mkdir -p ${dataDir}scplus_pipeline
scenicplus init_snakemake --out_dir ${dataDir}scplus_pipeline
mkdir -p ${dataDir}tmp

#### modify config file ####
vim ${dataDir}scplus_pipeline/Snakemake/config/config.yaml

## run Snakemake
# download genome annotation manually
cp -f /ix1/wchen/xiangyu/Ref_data/SCENIC/genome_annotation.tsv ${dataDir}scplus_pipeline/Snakemake/
cp -f /ix1/wchen/xiangyu/Ref_data/SCENIC/chromsizes.tsv ${dataDir}scplus_pipeline/Snakemake/
# Snakemake
cd ${dataDir}scplus_pipeline/Snakemake/
snakemake --cores 20

