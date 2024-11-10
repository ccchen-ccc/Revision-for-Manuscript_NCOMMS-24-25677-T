setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step2/data")
library(data.table)
library(parallel)
library(stringr)

mpra_vid=fread("sig05_case_AD_SC.vid",header=F)
names(mpra_vid)="vid"

#WGSanno#
anno=fread("sig05_case_AD_SC_WGSanno_extract")
mpra_vid=merge(mpra_vid,anno,all.x=TRUE,by="vid")

#reported cytoband#
report=fread("81_reported_loci")
mpra_vid$report_loci=ifelse(mpra_vid$cytoband %in% report$Cytoband,1,0)

#coreMarks#
chromHMM_active=c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts","10_TssBiv","11_BivFlnk","12_EnhBiv")
mpra_vid$chromHMM_active_tag=ifelse(mpra_vid$A549_coreMarks %in% chromHMM_active | mpra_vid$IMR90_coreMarks %in% chromHMM_active |mpra_vid$NHLF_coreMarks %in% chromHMM_active | mpra_vid$Fetallung_coreMarks %in% chromHMM_active |mpra_vid$Lung_coreMarks %in% chromHMM_active,1,0)

#coreMarks_epi#
mpra_vid$chromHMM_epi_active_tag=ifelse(mpra_vid$A549_coreMarks %in% chromHMM_active |mpra_vid$Fetallung_coreMarks %in% chromHMM_active |mpra_vid$Lung_coreMarks %in% chromHMM_active,1,0)

#eQTL_NMU#
eQTL_NMU=fread("sig05_case_AD_SC_vid_NMU_Lung_eqtl.txt")
mpra_vid=merge(mpra_vid,eQTL_NMU,by="vid",all.x=TRUE)

#添加EncodeLung注释#
EncodeLung=fread("sig05_case_AD_SC_vid_EncodeLung.bed",header=F)
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

#注释EncodeLung_epi_ATAC/Dnase#
peak_ATAC_Dnase_epi_tmp <- subset(EncodeLung,(grepl("ATAC",EncodeLung$V4)|grepl("DNase",EncodeLung$V4))& (grepl("lung",EncodeLung$V4)|grepl("PC-9",EncodeLung$V4)|grepl("epithelial",EncodeLung$V4)|grepl("H1299",EncodeLung$V4)|grepl("A549",EncodeLung$V4)|grepl("BEAS2B",EncodeLung$V4)))
mpra_vid$EncodeLung_epi_ATAC_Dnase_tag=ifelse(mpra_vid$vid %in% peak_ATAC_Dnase_epi_tmp$V8,1,0)

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

#注释EncodeLung_active histone2#
mpra_vid$EncodeLung_H3K27ac_H3K4me3_tag=ifelse(mpra_vid$EncodeLung_H3K27ac_tag==1 | mpra_vid$EncodeLung_H3K4me3_tag==1,1,0)

#注释EncodeLung_active histone1#
mpra_vid$EncodeLung_H3K27ac_H3K4me3_H3K4me1_H3K9ac_tag=ifelse(mpra_vid$EncodeLung_H3K27ac_tag==1 | mpra_vid$EncodeLung_H3K4me3_tag==1 | mpra_vid$EncodeLung_H3K4me1_tag==1 | mpra_vid$EncodeLung_H3K9ac_tag==1,1,0)

#注释EncodeLung_epi_H3K27ac#
peak_H3K27ac_tmp_epi <- subset(EncodeLung,(grepl("H3K27ac",EncodeLung$V4))& (grepl("lung",EncodeLung$V4)|grepl("PC-9",EncodeLung$V4)|grepl("epithelial",EncodeLung$V4)|grepl("H1299",EncodeLung$V4)|grepl("A549",EncodeLung$V4)|grepl("BEAS2B",EncodeLung$V4)))
mpra_vid$EncodeLung_epi_H3K27ac_tag=ifelse(mpra_vid$vid %in% peak_H3K27ac_tmp_epi$V8,1,0)

#注释EncodeLung_epi_H3K4me1#
peak_H3K4me1_tmp_epi <- subset(EncodeLung,(grepl("H3K4me1",EncodeLung$V4))& (grepl("lung",EncodeLung$V4)|grepl("PC-9",EncodeLung$V4)|grepl("epithelial",EncodeLung$V4)|grepl("H1299",EncodeLung$V4)|grepl("A549",EncodeLung$V4)|grepl("BEAS2B",EncodeLung$V4)))
mpra_vid$EncodeLung_epi_H3K4me1_tag=ifelse(mpra_vid$vid %in% peak_H3K4me1_tmp_epi$V8,1,0)

#注释EncodeLung_epi_H3K4me3#
peak_H3K4me3_tmp_epi <- subset(EncodeLung,(grepl("H3K4me3",EncodeLung$V4))& (grepl("lung",EncodeLung$V4)|grepl("PC-9",EncodeLung$V4)|grepl("epithelial",EncodeLung$V4)|grepl("H1299",EncodeLung$V4)|grepl("A549",EncodeLung$V4)|grepl("BEAS2B",EncodeLung$V4)))
mpra_vid$EncodeLung_epi_H3K4me3_tag=ifelse(mpra_vid$vid %in% peak_H3K4me3_tmp_epi$V8,1,0)

#注释EncodeLung_epi_H3K9ac#
peak_H3K9ac_tmp_epi <- subset(EncodeLung,(grepl("H3K9ac",EncodeLung$V4))& (grepl("lung",EncodeLung$V4)|grepl("PC-9",EncodeLung$V4)|grepl("epithelial",EncodeLung$V4)|grepl("H1299",EncodeLung$V4)|grepl("A549",EncodeLung$V4)|grepl("BEAS2B",EncodeLung$V4)))
mpra_vid$EncodeLung_epi_H3K9ac_tag=ifelse(mpra_vid$vid %in% peak_H3K9ac_tmp_epi$V8,1,0)

#注释EncodeLung_active histone2#
mpra_vid$EncodeLung_epi_H3K27ac_H3K4me3_tag=ifelse(mpra_vid$EncodeLung_epi_H3K27ac_tag==1 | mpra_vid$EncodeLung_epi_H3K4me3_tag==1,1,0)

#注释EncodeLung_active histone1#
mpra_vid$EncodeLung_epi_H3K27ac_H3K4me3_H3K4me1_H3K9ac_tag=ifelse(mpra_vid$EncodeLung_epi_H3K27ac_tag==1 | mpra_vid$EncodeLung_epi_H3K4me3_tag==1 | mpra_vid$EncodeLung_epi_H3K4me1_tag==1 | mpra_vid$EncodeLung_epi_H3K9ac_tag==1,1,0)

#注释EncodeLung_H3K36me3#
peak_H3K36me3_tmp <- subset(EncodeLung,grepl("H3K36me3",EncodeLung$V4))
mpra_vid$EncodeLung_H3K36me3_tag=ifelse(mpra_vid$vid %in% peak_H3K36me3_tmp$V8,1,0)

#注释EncodeLung_H3K9me3#
peak_H3K9me3_tmp <- subset(EncodeLung,grepl("H3K9me3",EncodeLung$V4))
mpra_vid$EncodeLung_H3K9me3_tag=ifelse(mpra_vid$vid %in% peak_H3K9me3_tmp$V8,1,0)

#注释EncodeLung_H3K27me3#
peak_H3K27me3_tmp <- subset(EncodeLung,grepl("H3K27me3",EncodeLung$V4))
mpra_vid$EncodeLung_H3K27me3_tag=ifelse(mpra_vid$vid %in% peak_H3K27me3_tmp$V8,1,0)

#注释EncodeLung_inactive histone1#
mpra_vid$EncodeLung_H3K27me3_H3K36me3_H3K9me3_tag=ifelse(mpra_vid$EncodeLung_H3K27me3_tag==1 | mpra_vid$EncodeLung_H3K36me3_tag==1 | mpra_vid$EncodeLung_H3K9me3_tag==1,1,0)

#注释EncodeLung_CTCF#
peak_CTCF_tmp <- subset(EncodeLung,grepl("CTCF",EncodeLung$V4))
mpra_vid$EncodeLung_CTCF_tag=ifelse(mpra_vid$vid %in% peak_CTCF_tmp$V8,1,0)

#注释EncodeLung_EP300#
peak_EP300_tmp <- subset(EncodeLung,grepl("EP300",EncodeLung$V4))
mpra_vid$EncodeLung_EP300_tag=ifelse(mpra_vid$vid %in% peak_EP300_tmp$V8,1,0)

#注释EncodeLung_POLR2A#
peak_POLR2A_tmp <- subset(EncodeLung,grepl("POLR2A",EncodeLung$V4))
mpra_vid$EncodeLung_POLR2A_tag=ifelse(mpra_vid$vid %in% peak_POLR2A_tmp$V8,1,0)


#注释LungENN#
deephaem=fread("./LungENN/prediction/sig05_case_AD_SC_vid_LungENN.txt")
mpra_vid=merge(mpra_vid,deephaem,all.x=TRUE,by.x="vid",by.y="name")

#注释TFBS_motifbreakR#
TFBS_motifbreakR=fread("TFBS_motifbreakR.txt")
mpra_vid=merge(mpra_vid,TFBS_motifbreakR,all.x=TRUE,by="vid")

#注释TFBS_SNP2TFBS#
SNP2TFBS=fread("TFBS_SNP2TFBS.txt")
SNP2TFBS=subset(SNP2TFBS,ScoreDiff!=0)[,c("vid","TFBS_SNP2TFBS")]
SNP2TFBS=SNP2TFBS[order(SNP2TFBS$vid),]

TF=unlist(sapply(1:length(SNP2TFBS$vid),function(x){str_c(subset(SNP2TFBS,vid==SNP2TFBS$vid[x])$TFBS_SNP2TFBS,collapse=";")}))
vid=SNP2TFBS$vid

TFBS_SNP2TFBS=unique(data.frame("vid"=vid,"TFBS_SNP2TFBS"=TF))
TFBS_SNP2TFBS=as.data.frame(TFBS_SNP2TFBS)

mpra_vid=merge(mpra_vid,TFBS_SNP2TFBS,all.x=TRUE,by="vid")

write.table(mpra_vid,"sig05_case_AD_SC_vid.anno",row.names=F,col.names=T,sep="\t",quote=F)
