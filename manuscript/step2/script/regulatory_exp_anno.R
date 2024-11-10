setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step2/data")
library(data.table)
library(parallel)
library(stringr)

mpra_vid=fread("sig05_case_AD_SC.vid",header=F)
names(mpra_vid)="vid"

#添加EncodeLung注释#
EncodeLung=fread("sig05_case_AD_SC_60bp_EncodeLung.bed",header=F)
EncodeLung=EncodeLung[,c("V8","V4")]

peak=unlist(sapply(1:length(EncodeLung$V8),function(x){str_c(subset(EncodeLung,V8==EncodeLung$V8[x])$V4,collapse=";")}))
vid=EncodeLung$V8
EncodeLung_peak=unique(data.frame("vid"=vid,"EncodeLung_peak"=peak))
EncodeLung_peak=as.data.frame(EncodeLung_peak)

mpra_vid=merge(mpra_vid,EncodeLung_peak,by="vid",all.x=TRUE)

#注释EncodeLung_DNase#
peak_DNase_tmp <- subset(EncodeLung,grepl("DNase",EncodeLung$V4))
mpra_vid$EncodeLung_DNase_tag=ifelse(mpra_vid$vid %in% peak_DNase_tmp$V8,1,0)

#注释EncodeLung_ATAC#
peak_ATAC_tmp <- subset(EncodeLung,grepl("ATAC",EncodeLung$V4))
mpra_vid$EncodeLung_ATAC_tag=ifelse(mpra_vid$vid %in% peak_ATAC_tmp$V8,1,0)

#注释EncodeLung_ATAC/Dnase#
peak_ATAC_Dnase_tmp <- subset(EncodeLung,grepl("ATAC",EncodeLung$V4)|grepl("DNase",EncodeLung$V4))
mpra_vid$EncodeLung_ATAC_Dnase_tag=ifelse(mpra_vid$vid %in% peak_ATAC_Dnase_tmp$V8,1,0)

#注释EncodeLung_H3K27ac#
peak_H3K27ac_tmp <- subset(EncodeLung,grepl("H3K27ac",EncodeLung$V4))
mpra_vid$EncodeLung_H3K27ac_tag=ifelse(mpra_vid$vid %in% peak_H3K27ac_tmp$V8,1,0)

#注释EncodeLung_H3K4me3#
peak_H3K4me3_tmp <- subset(EncodeLung,grepl("H3K4me3",EncodeLung$V4))
mpra_vid$EncodeLung_H3K4me3_tag=ifelse(mpra_vid$vid %in% peak_H3K4me3_tmp$V8,1,0)

#注释EncodeLung_H3K4me1#
peak_H3K4me1_tmp <- subset(EncodeLung,grepl("H3K4me1",EncodeLung$V4))
mpra_vid$EncodeLung_H3K4me1_tag=ifelse(mpra_vid$vid %in% peak_H3K4me1_tmp$V8,1,0)

#注释EncodeLung_H3K9ac#
peak_H3K9ac_tmp <- subset(EncodeLung,grepl("H3K9ac",EncodeLung$V4))
mpra_vid$EncodeLung_H3K9ac_tag=ifelse(mpra_vid$vid %in% peak_H3K9ac_tmp$V8,1,0)

#注释EncodeLung_H3K36me3#
peak_H3K36me3_tmp <- subset(EncodeLung,grepl("H3K36me3",EncodeLung$V4))
mpra_vid$EncodeLung_H3K36me3_tag=ifelse(mpra_vid$vid %in% peak_H3K36me3_tmp$V8,1,0)

#注释EncodeLung_H3K9me3#
peak_H3K9me3_tmp <- subset(EncodeLung,grepl("H3K9me3",EncodeLung$V4))
mpra_vid$EncodeLung_H3K9me3_tag=ifelse(mpra_vid$vid %in% peak_H3K9me3_tmp$V8,1,0)

#注释EncodeLung_H3K27ac_clean#
mpra_vid$EncodeLung_H3K27ac_clean_tag=ifelse(mpra_vid$EncodeLung_H3K27ac_tag==1 & mpra_vid$EncodeLung_H3K36me3_tag!=1 & mpra_vid$EncodeLung_H3K9me3_tag!=1,1,0)

#注释EncodeLung_H3K4me3_clean#
mpra_vid$EncodeLung_H3K4me3_clean_tag=ifelse(mpra_vid$EncodeLung_H3K4me3_tag==1 & mpra_vid$EncodeLung_H3K36me3_tag!=1 & mpra_vid$EncodeLung_H3K9me3_tag!=1,1,0)

#注释注释EncodeLung_H3K36me3_clean#
mpra_vid$EncodeLung_H3K36me3_clean_tag=ifelse(mpra_vid$EncodeLung_H3K36me3_tag==1 & mpra_vid$EncodeLung_H3K27ac_tag!=1 & mpra_vid$EncodeLung_H3K4me3_tag!=1,1,0)

#注释注释EncodeLung_H3K9me3_clean#
mpra_vid$EncodeLung_H3K9me3_clean_tag=ifelse(mpra_vid$EncodeLung_H3K9me3_tag==1 & mpra_vid$EncodeLung_H3K27ac_tag!=1 & mpra_vid$EncodeLung_H3K4me3_tag!=1,1,0)

#注释EncodeLung_CTCF#
peak_CTCF_tmp <- subset(EncodeLung,grepl("CTCF",EncodeLung$V4))
mpra_vid$EncodeLung_CTCF_tag=ifelse(mpra_vid$vid %in% peak_CTCF_tmp$V8,1,0)

#注释EncodeLung_EP300#
peak_EP300_tmp <- subset(EncodeLung,grepl("EP300",EncodeLung$V4))
mpra_vid$EncodeLung_EP300_tag=ifelse(mpra_vid$vid %in% peak_EP300_tmp$V8,1,0)

#注释EncodeLung_POLR2A#
peak_POLR2A_tmp <- subset(EncodeLung,grepl("POLR2A",EncodeLung$V4))
mpra_vid$EncodeLung_POLR2A_tag=ifelse(mpra_vid$vid %in% peak_POLR2A_tmp$V8,1,0)

#注释EncodeLung_A549_DNase#
peak_A549_DNase_tmp <- subset(EncodeLung,grepl("A549_DNase",EncodeLung$V4))
mpra_vid$EncodeLung_A549_DNase_tag=ifelse(mpra_vid$vid %in% peak_A549_DNase_tmp$V8,1,0)

#注释EncodeLung_A549_ATAC#
peak_A549_ATAC_tmp <- subset(EncodeLung,grepl("A549_ATAC",EncodeLung$V4))
mpra_vid$EncodeLung_A549_ATAC_tag=ifelse(mpra_vid$vid %in% peak_A549_ATAC_tmp$V8,1,0)

#注释EncodeLung_A549_ATAC/Dnase#
peak_A549_ATAC_Dnase_tmp <- subset(EncodeLung,grepl("A549_ATAC",EncodeLung$V4)|grepl("A549_DNase",EncodeLung$V4))
mpra_vid$EncodeLung_A549_ATAC_Dnase_tag=ifelse(mpra_vid$vid %in% peak_A549_ATAC_Dnase_tmp$V8,1,0)

#注释EncodeLung_A549_H3K27ac#
peak_A549_H3K27ac_tmp <- subset(EncodeLung,grepl("A549_H3K27ac",EncodeLung$V4))
mpra_vid$EncodeLung_A549_H3K27ac_tag=ifelse(mpra_vid$vid %in% peak_A549_H3K27ac_tmp$V8,1,0)

#注释EncodeLung_A549_H3K4me3#
peak_A549_H3K4me3_tmp <- subset(EncodeLung,grepl("A549_H3K4me3",EncodeLung$V4))
mpra_vid$EncodeLung_A549_H3K4me3_tag=ifelse(mpra_vid$vid %in% peak_A549_H3K4me3_tmp$V8,1,0)

#注释EncodeLung_A549_H3K4me1#
peak_A549_H3K4me1_tmp <- subset(EncodeLung,grepl("A549_H3K4me1",EncodeLung$V4))
mpra_vid$EncodeLung_A549_H3K4me1_tag=ifelse(mpra_vid$vid %in% peak_A549_H3K4me1_tmp$V8,1,0)

#注释EncodeLung_A549_H3K9ac#
peak_A549_H3K9ac_tmp <- subset(EncodeLung,grepl("A549_H3K9ac",EncodeLung$V4))
mpra_vid$EncodeLung_A549_H3K9ac_tag=ifelse(mpra_vid$vid %in% peak_A549_H3K9ac_tmp$V8,1,0)

#注释EncodeLung_A549_H3K36me3#
peak_A549_H3K36me3_tmp <- subset(EncodeLung,grepl("A549_H3K36me3",EncodeLung$V4))
mpra_vid$EncodeLung_A549_H3K36me3_tag=ifelse(mpra_vid$vid %in% peak_A549_H3K36me3_tmp$V8,1,0)

#注释EncodeLung_A549_H3K9me3#
peak_A549_H3K9me3_tmp <- subset(EncodeLung,grepl("A549_H3K9me3",EncodeLung$V4))
mpra_vid$EncodeLung_A549_H3K9me3_tag=ifelse(mpra_vid$vid %in% peak_A549_H3K9me3_tmp$V8,1,0)

#注释EncodeLung_A549_CTCF#
peak_A549_CTCF_tmp <- subset(EncodeLung,grepl("A549_CTCF",EncodeLung$V4))
mpra_vid$EncodeLung_A549_CTCF_tag=ifelse(mpra_vid$vid %in% peak_A549_CTCF_tmp$V8,1,0)

#注释EncodeLung_A549_EP300#
peak_A549_EP300_tmp <- subset(EncodeLung,grepl("A549_EP300",EncodeLung$V4))
mpra_vid$EncodeLung_A549_EP300_tag=ifelse(mpra_vid$vid %in% peak_A549_EP300_tmp$V8,1,0)

#注释EncodeLung_A549_POLR2A#
peak_A549_POLR2A_tmp <- subset(EncodeLung,grepl("A549_POLR2A",EncodeLung$V4))
mpra_vid$EncodeLung_A549_POLR2A_tag=ifelse(mpra_vid$vid %in% peak_A549_POLR2A_tmp$V8,1,0)

#注释EncodeLung_H1299_DNase#
peak_H1299_DNase_tmp <- subset(EncodeLung,grepl("H1299_DNase",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_DNase_tag=ifelse(mpra_vid$vid %in% peak_H1299_DNase_tmp$V8,1,0)

#注释EncodeLung_H1299_ATAC#
peak_H1299_ATAC_tmp <- subset(EncodeLung,grepl("H1299_ATAC",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_ATAC_tag=ifelse(mpra_vid$vid %in% peak_H1299_ATAC_tmp$V8,1,0)

#注释EncodeLung_H1299_ATAC/Dnase#
peak_H1299_ATAC_Dnase_tmp <- subset(EncodeLung,grepl("H1299_ATAC",EncodeLung$V4)|grepl("H1299_DNase",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_ATAC_Dnase_tag=ifelse(mpra_vid$vid %in% peak_H1299_ATAC_Dnase_tmp$V8,1,0)

#注释EncodeLung_H1299_H3K27ac#
peak_H1299_H3K27ac_tmp <- subset(EncodeLung,grepl("H1299_H3K27ac",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_H3K27ac_tag=ifelse(mpra_vid$vid %in% peak_H1299_H3K27ac_tmp$V8,1,0)

#注释EncodeLung_H1299_H3K4me3#
peak_H1299_H3K4me3_tmp <- subset(EncodeLung,grepl("H1299_H3K4me3",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_H3K4me3_tag=ifelse(mpra_vid$vid %in% peak_H1299_H3K4me3_tmp$V8,1,0)

#注释EncodeLung_H1299_H3K4me1#
peak_H1299_H3K4me1_tmp <- subset(EncodeLung,grepl("H1299_H3K4me1",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_H3K4me1_tag=ifelse(mpra_vid$vid %in% peak_H1299_H3K4me1_tmp$V8,1,0)

#注释EncodeLung_H1299_H3K9ac#
peak_H1299_H3K9ac_tmp <- subset(EncodeLung,grepl("H1299_H3K9ac",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_H3K9ac_tag=ifelse(mpra_vid$vid %in% peak_H1299_H3K9ac_tmp$V8,1,0)

#注释EncodeLung_H1299_H3K36me3#
peak_H1299_H3K36me3_tmp <- subset(EncodeLung,grepl("H1299_H3K36me3",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_H3K36me3_tag=ifelse(mpra_vid$vid %in% peak_H1299_H3K36me3_tmp$V8,1,0)

#注释EncodeLung_H1299_H3K9me3#
peak_H1299_H3K9me3_tmp <- subset(EncodeLung,grepl("H1299_H3K9me3",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_H3K9me3_tag=ifelse(mpra_vid$vid %in% peak_H1299_H3K9me3_tmp$V8,1,0)

#注释EncodeLung_H1299_CTCF#
peak_H1299_CTCF_tmp <- subset(EncodeLung,grepl("H1299_CTCF",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_CTCF_tag=ifelse(mpra_vid$vid %in% peak_H1299_CTCF_tmp$V8,1,0)

#注释EncodeLung_H1299_EP300#
peak_H1299_EP300_tmp <- subset(EncodeLung,grepl("H1299_EP300",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_EP300_tag=ifelse(mpra_vid$vid %in% peak_H1299_EP300_tmp$V8,1,0)

#注释EncodeLung_H1299_POLR2A#
peak_H1299_POLR2A_tmp <- subset(EncodeLung,grepl("H1299_POLR2A",EncodeLung$V4))
mpra_vid$EncodeLung_H1299_POLR2A_tag=ifelse(mpra_vid$vid %in% peak_H1299_POLR2A_tmp$V8,1,0)

write.table(mpra_vid,"sig05_case_AD_SC_60bp.anno",row.names=F,col.names=T,sep="\t",quote=F)
