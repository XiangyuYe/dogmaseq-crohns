#!/bin/bash
#
#SBATCH --time=12:00:00
#SBATCH --job-name=dogma_qc
#SBATCH --mem=100G
#SBATCH --cpus-per-task=10

#SBATCH --array=1-42
#SBATCH --output=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/00_temp/00_dogma_qc%a.out
#SBATCH --error=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/00_temp/00_dogma_qc%a.err

bash
let k=0

## data path
proj_path=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/
clean_path=${proj_path}02_clean/
cr_path=/ix1/rduerr/shared/rduerr_wchen/raw/Cell_Ranger_output/cellranger-arc_count_2.0.2/
feature_ref=/ix1/rduerr/shared/rduerr_wchen/Experiment_summary_feature_ref_cite_seq_index/feature_ref_TotalSeq-A-Human-Universal-Cocktail-V1_HTOs.csv
kite_path=/ix1/wchen/xiangyu/Biosoft/kite/
sample_list=${proj_path}Helmsley_sample.txt
n_sample=`cat ${sample_list} | wc -l`

for fdx in `seq 1 ${n_sample}`
do

let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
then

idx=`head -n ${fdx} ${sample_list} | tail -n 1 | awk -F"\t" '{print $1}'`
batch_fqx=`head -n ${fdx} ${sample_list} | tail -n 1 | awk -F"\t" '{print $2}'`
batch_crx=`head -n ${fdx} ${sample_list} | tail -n 1 | awk -F"\t" '{print $3}'`

#
in_path_arc=${cr_path}${batch_crx}/${idx}_arc/
in_path_adt=${proj_path}01_raw/${idx}/ADT/featurecounts/
in_path_hto=${proj_path}01_raw/${idx}/HTO/featurecounts/
hto_valid_file=${proj_path}01_raw/${idx}/HTO/Feature_HTOs_valid.csv

# # renames cells in fragment file
source activate /ix1/wchen/xiangyu/conda_env/scanpy

frag_path=/ix1/rduerr/shared/rduerr_wchen/Xiangyu/Helmsley/Fragment_modified/${idx}/
mkdir -p ${frag_path}
zcat ${in_path_arc}atac_fragments.tsv.gz | 
awk -v prefix="$idx" \
'BEGIN {OFS="\t"} /^#/ {print $0} !/^#/ {print $1, $2, $3, prefix "_" $4, $5}' | \
bgzip > ${frag_path}modified_fragments.tsv.gz
tabix -p bed ${frag_path}modified_fragments.tsv.gz

conda deactivate

# qun qc
qc_R=${proj_path}code/00_01_04_QC.R
module load r/4.5.0
Rscript ${qc_R} \
--in_path_arc ${in_path_arc} \
--arc_fragfile ${proj_path}01_raw/${idx}/modified_fragments.tsv.gz \
--in_path_adt ${in_path_adt} \
--in_path_hto ${in_path_hto} \
--hto_valid_file ${hto_valid_file} \
--clean_path ${clean_path} \
--threads 10 \
--idx ${idx}

module unload r/4.5.0

fi
done

