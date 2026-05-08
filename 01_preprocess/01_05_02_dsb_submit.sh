#!/bin/bash
#
#SBATCH --time=12:00:00
#SBATCH --job-name=dogma_dsb
#SBATCH --mem=100G
#SBATCH --cpus-per-task=10

#SBATCH --array=1-42
#SBATCH --output=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/00_temp/00_dogma_dsb%a.out
#SBATCH --error=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/00_temp/00_dogma_dsb%a.err

bash
let k=0

## data path
proj_path=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/
clean_path=${proj_path}02_clean/
cr_path=/ix1/rduerr/shared/rduerr_wchen/raw/Cell_Ranger_output/cellranger-arc_count_2.0.2/
dsb_list=${proj_path}Helmsley_pre_dsb.txt
n_sample=`cat ${dsb_list} | wc -l`

for fdx in `seq 1 ${n_sample}`
do

let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
then

idx=`head -n ${fdx} ${dsb_list} | tail -n 1 | awk -F"\t" '{print $1}'`
batch_crx=`head -n ${fdx} ${dsb_list} | tail -n 1 | awk -F"\t" '{print $3}'`
low_prot=`head -n ${fdx} ${dsb_list} | tail -n 1 | awk -F"\t" '{print $4}'`
hi_prot=`head -n ${fdx} ${dsb_list} | tail -n 1 | awk -F"\t" '{print $5}'`
#
in_path_arc=${cr_path}${batch_crx}/${idx}_arc/
in_path_adt=${proj_path}01_raw/${idx}/ADT/featurecounts/

# qun qc
dsb_adt=${proj_path}code/00_01_05_dsb_adt.R
module load r/4.5.0
Rscript ${dsb_adt} \
--in_path_arc ${in_path_arc} \
--in_path_adt ${in_path_adt} \
--clean_path ${clean_path} \
--low_prot ${low_prot} \
--hi_prot ${hi_prot} \
--idx ${idx}

module unload r/4.5.0

fi
done

