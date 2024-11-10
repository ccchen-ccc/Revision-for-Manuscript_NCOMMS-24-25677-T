setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step2/data")
library(data.table)
library(parallel)
library(stringr)
args = c("3cell_cov_rep6",0.01)

#MPRA#
mpra=fread("vid_alleleEffect_for_MPRA.txt")

#anno#
anno=fread("sig05_case_AD_SC_vid.anno")
mpra=merge(mpra,anno,by="vid",all.x=TRUE)

#GWAS_asso#
asso=fread("sig05_case_AD_SC.asso")

mpra=merge(mpra,asso,all.X=TRUE,by="vid")

write.table(mpra,paste("sig05_case_AD_SC_vid_alleleEffect_",args[1],"_",as.numeric(args[2]),"_filter.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
