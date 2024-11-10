setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step2/data")
library(data.table)
library(parallel)
library(stringr)
args = c("3cell_cov_rep6",0.01)

#MPRA效应#
mpra=fread("vid_alleleEffect_for_MPRA.txt")

#anno#
anno=fread("sig05_case_AD_SC_60bp.anno")
mpra=merge(mpra,anno,by="vid",all.x=TRUE)

write.table(mpra,"sig05_case_AD_SC_60bp_expEffect_3cell_cov_rep6_0.01_filter.txt",row.names=F,col.names=T,sep="\t",quote=F)
