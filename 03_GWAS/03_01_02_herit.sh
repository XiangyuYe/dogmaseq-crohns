
ldsc=/ix1/wchen/xiangyu/Biosoft/ldsc/ldsc.py
summ_path=/ix1/wchen/xiangyu/Ref_data/GWAS/GRCh37/
out_path=/ix1/wchen/xiangyu/Ref_data/GWAS/GRCh37/h2/
ref_path=/ix1/wchen/xiangyu/Ref_data/ldsc/

p=CD
sum_prefix=${summ_path}${p}_hm3
out_prefix=${out_path}${p}_hm3

python ${ldsc} \
--h2 ${sum_prefix}.ldsc.gz \
--ref-ld-chr ${ref_path}eur_ref_ld_chr/ \
--w-ld-chr ${ref_path}eur_w_ld_chr/ \
--out ${out_prefix}_h2

