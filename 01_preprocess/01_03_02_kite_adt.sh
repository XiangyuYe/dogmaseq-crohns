#!/bin/bash
#
#SBATCH --time=12:00:00
#SBATCH --job-name=dogma_qc1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=20

#SBATCH --array=1-42
#SBATCH --output=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/00_temp/00_dogma_adt%a.out
#SBATCH --error=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/00_temp/00_dogma_adt%a.err

bash
let k=0

## data path
proj_path=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/
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

echo ${idx}

raw_path=/ix1/rduerr/shared/rduerr_wchen/raw/Single_cell_fastq/${batch_fqx}_fastq/${idx}_fastq/
out_path=${proj_path}01_raw/${idx}/ADT/
mkdir -p ${out_path}
##
arc_out_path=${cr_path}${batch_crx}${batchx}/${idx}_arc/
cut -f 1 -d , ${arc_out_path}/per_barcode_metrics.csv | \
sed 's/-1//g' | grep barcode -v > ${out_path}barcodes.txt
#
cat ${feature_ref} | grep "ADT_" | \
awk -F',' 'BEGIN {OFS=","} {print $2,$5}' > ${out_path}Feature_ADTs.csv
##
module load python/3.7.0 kallisto/0.50.1 bustools/0.45.1
python ${kite_path}featuremap/featuremap.py \
${out_path}Feature_ADTs.csv \
--t2g ${out_path}FeaturesMismatch.t2g \
--fa ${out_path}FeaturesMismatch.fa 
#
kallisto index \
-i ${out_path}FeaturesMismatch.idx \
-k 15 \
${out_path}FeaturesMismatch.fa
#
kallisto bus \
-i ${out_path}FeaturesMismatch.idx \
-o ${out_path} \
-x 10xv3 \
-t 24 \
${raw_path}*ADT*
#
bustools correct \
-w ${out_path}barcodes.txt \
${out_path}output.bus \
-o ${out_path}output_corrected.bus
#
bustools sort \
-t 4 \
-o ${out_path}output_sorted.bus ${out_path}output_corrected.bus
#
mkdir -p ${out_path}featurecounts/
bustools count \
-o ${out_path}featurecounts/featurecounts \
--genecounts \
-g ${out_path}FeaturesMismatch.t2g \
-e ${out_path}matrix.ec \
-t ${out_path}transcripts.txt ${out_path}output_sorted.bus

module unload python/3.7.0 kallisto/0.50.1 bustools/0.45.1

fi
done
