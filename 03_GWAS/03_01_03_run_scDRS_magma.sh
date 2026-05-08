## 0. set parameters
grch=37.3
eth=eur
#
project_path=/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/
out_path=${project_path}03_output/01_scDRS/
#
magma=/ix1/wchen/xiangyu/Biosoft/MAGMA/magma
zscore_format=${project_path}code/FUNCTION/scDRS/00_scDRS_format_zscore.R
ref_path=/ix1/wchen/xiangyu/Ref_data/

mkdir -p ${out_path}/01_magma
mkdir -p ${out_path}/02_gs
mkdir -p ${out_path}/03_integrate
mkdir -p ${out_path}/04_downstream

## 1. make gene annotation file
conda activate /ix1/wchen/xiangyu/conda_env/scanpy
module load r/4.5.0
${magma} \
--annotate window=10,10 \
--snp-loc ${ref_path}genome/homo_sapiens/g1000_eur.bim \
--gene-loc ${ref_path}genome/homo_sapiens/NCBI${grch}.gene.loc \
--out ${out_path}/01_magma/step1

# 2. run MAGMA and scDRS for each summary (may submit to PBS)
# 2.1 set summary
for p in CD IBD UC
do
n_top=1000
n_ctrl=1000
summ_file=${ref_path}GWAS/GRCh37/${p}.txt
sc_path=${project_path}/03_output/01_scDRS/00_scRNA/
## 2.2 magma analysis
${magma} \
--bfile ${ref_path}/genome/homo_sapiens/g1000_${eth} \
--pval ${summ_file} use='rsid,p' ncol='n' \
--gene-annot ${out_path}/01_magma/step1.genes.annot \
--out ${out_path}/01_magma/${p}
#
Rscript ${zscore_format} \
--megama_out ${out_path}/01_magma/${p}.genes.out \
--grch 37 \
--pheno ${p}

## 2.3 scDRS analysis
# Select top 1,000 genes and use z-score weights
scdrs munge-gs \
--out-file ${out_path}/02_gs/${p}_top${n_top}.gs \
--zscore-file ${out_path}/01_magma/zscore_${p}.tsv \
--weight zscore \
--n-max ${n_top}

done

conda deactivate
