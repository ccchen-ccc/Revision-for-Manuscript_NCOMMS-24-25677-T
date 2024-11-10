library(data.table)
library(stringr)
setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step2/data")

vid=fread("sig05_case_AD_SC.vid",header=F)
anno=fread("genecodev19_anno.txt")
anno=subset(anno,type=="protein_coding")
anno=anno[,c("gene_id","symbol","ENSG")]

#NMU_Lung0.05#
eqtl=fread("/public/home/ccc2018/data/database/LC_WGS_germline/eqtl_v2/normal/result/nominal/lung_normal_eqtl.sigpairs.txt")
dat=subset(eqtl,variant_id %in% vid$V1)
dat$ref_mpra=sapply(strsplit(dat$variant_id,':'),'[',3)
dat$slope_adj=ifelse(dat$ref_mpra==dat$ref,dat$slope,-1*dat$slope)
dat1=dat[,c("variant_id","gene_id","slope_adj")]

dat2=merge(dat1,anno)

eqtl_gene=unlist(sapply(1:length(dat2$variant_id),function(x){str_c(subset(dat2,variant_id==dat2$variant_id[x])$symbol,collapse=";")}))
eqtl_gene_slope=unlist(sapply(1:length(dat2$variant_id),function(x){str_c(round(subset(dat2,variant_id==dat2$variant_id[x])$slope_adj,3),collapse=";")}))

vid_eqtl=dat2$variant_id
dat3=unique(data.frame("vid"=vid_eqtl,"NMU_Lung_eqtl_gene"=eqtl_gene,"NMU_Lung_eqtl_gene_slope"=eqtl_gene_slope))


write.table(dat3,"sig05_case_AD_SC_vid_NMU_Lung_eqtl.txt",row.names=F,col.names=T,sep="\t",quote=F)




