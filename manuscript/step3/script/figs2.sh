#############Example: findMotifsGenome.pl peaks.txt mm8r peakAnalysis -size 200 -len 8
###
tools=/public/home/xj2021/software/homer/bin/findMotifsGenome.pl
expr_pos=/public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs2/oligo_expr_pos.txt
backgroud=/public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs2/oligo_expr_background.txt
pos_out=./pos_homer
expr_neg=/public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs2/oligo_expr_neg.txt
neg_out=./neg_homer
perl ${tools} ${expr_pos} hg19 ${pos_out} -bg ${backgroud}& 
perl ${tools} ${expr_neg} hg19 ${neg_out} -bg ${backgroud}&