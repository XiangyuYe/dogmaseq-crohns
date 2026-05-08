#!/bin/bash
#
#SBATCH --time=3-00:00:00
#SBATCH --job-name=cisTarget
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20

#SBATCH --array=1
#SBATCH --output=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/00_temp/cisTarget%a.out
#SBATCH --error=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/00_temp/cisTarget%a.err

bash

## input parameters
outDir=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/03_output/07_SCENIC/
region_bed=${outDir}consensus_regions.bed
out_prefix="CD"
outfa=${outDir}hg38.${out_prefix}.with_1kb_bg_padding.fa

## fixed parameters
ref_path=/ix1/wchen/xiangyu/Ref_data/
soft_dir=/ix1/wchen/xiangyu/Biosoft/create_cisTarget_databases/

## set data file path
genome_fasta=${ref_path}genome/homo_sapiens/UCSC_hg38/hg38.fa
chromsizes=${ref_path}genome/homo_sapiens/UCSC_hg38/hg38.chrom.sizes
cb_dir=${ref_path}SCENIC/aertslab_motif_colleciton/v10nr_clust_public/singletons/
ls ${cb_dir} > ${outDir}motifs.txt

source activate /ix1/wchen/xiangyu/conda_env/scenicplus
## format bed file
${soft_dir}create_fasta_with_padded_bg_from_bed.sh \
${genome_fasta} \
${chromsizes} \
${region_bed} \
${outfa} \
1000 \
yes

## run cistarget
${soft_dir}/create_cistarget_motif_databases.py \
-f ${outfa} \
-M ${cb_dir} \
-m ${outDir}motifs.txt \
-o ${outDir}/${out_prefix} \
--bgpadding 1000 \
-t 20

conda deactivate

