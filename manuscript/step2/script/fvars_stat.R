library(data.table)
library(parallel)

setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step2/data")

dat=fread("sig05_case_AD_SC_vid_alleleEffect_3cell_cov_rep6_0.01_filter.txt")
dat1=dat[order(dat$vid,-dat$padj_allele_tag,-dat$padj_scram_tag,dat$padj_allele),]
dat1=dat1[!duplicated(dat1$vid),]

dat1$Type=ifelse(dat1$padj_allele_tag==1 & dat1$padj_scram_tag==1,"Fvars","Non-fvars")

stat=as.data.frame(table(dat1$cytoband))

dat1_2=unique(subset(dat1,Type=="Fvars" )[,c("vid","Type","ref","alt","rsID","chr","pos","region","cytoband","AF_CN_LC","P_LC","log2FoldChange_allele","report_loci","A549_coreMarks","IMR90_coreMarks","NHLF_coreMarks","Fetallung_coreMarks","Lung_coreMarks","chromHMM_active_tag","EncodeLung_peak","EncodeLung_ATAC_Dnase_tag","EncodeLung_H3K27ac_H3K4me3_H3K4me1_H3K9ac_tag","NMU_Lung_eqtl_gene","TFBS_HOCOMOCO","TFBS_jASPAR","TFBS_SNP2TFBS","EncodeLung_ATAC_Dnase_tag","EncodeLung_H3K27ac_H3K4me3_H3K4me1_H3K9ac_tag")])

#fvars region stat#
func=function(i){
print (i)
res=data.frame(Chr=subset(dat1,cytoband==i)$chr[1],locus=i,Lead_snp=subset(dat1,P_LC==min(subset(dat1,cytoband==i)$P_LC))$vid,No_of_SNPs_tested=nrow(subset(dat1,cytoband==i)),No_of_Fvars=nrow(subset(dat1_2,cytoband==i & Type=="Fvars")))
return(res)
}
i<- unique(dat1_2$cytoband)
cl <- makeCluster(10,type="FORK")
num <- parLapply(cl,i,func)
library(plyr)
results<-ldply(num,data.frame)
stopCluster(cl)
results=results[order(results$Chr),]

write.table(results,"fvars_region_stat.txt",row.names=F,col.names=T,quote=F,sep="\t")

dat2=unique(subset(dat1,Type=="Fvars" & (EncodeLung_ATAC_Dnase_tag==1 | EncodeLung_H3K27ac_H3K4me3_H3K4me1_H3K9ac_tag==1) & !((is.na(NMU_Lung_eqtl_gene) & is.na(TFBS_HOCOMOCO) & is.na(TFBS_jASPAR) & is.na(TFBS_SNP2TFBS))))[,c("vid","Type","ref","alt","rsID","chr","pos","region","cytoband","AF_CN_LC","P_LC","log2FoldChange_allele","report_loci","A549_coreMarks","IMR90_coreMarks","NHLF_coreMarks","Fetallung_coreMarks","Lung_coreMarks","chromHMM_active_tag","EncodeLung_peak","EncodeLung_ATAC_Dnase_tag","EncodeLung_H3K27ac_H3K4me3_H3K4me1_H3K9ac_tag","NMU_Lung_eqtl_gene","TFBS_HOCOMOCO","TFBS_jASPAR","TFBS_SNP2TFBS","EncodeLung_ATAC_Dnase_tag","EncodeLung_H3K27ac_H3K4me3_H3K4me1_H3K9ac_tag")])

#mpra_vars#
for (i in unique(dat2$cytoband)){
tmp=subset(dat1,cytoband==i)[,c("vid")]
write.table(tmp,paste(i,".mpravars",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
}

#fvars region stat#
func=function(i){
print (i)
res=data.frame(Chr=subset(dat1,cytoband==i)$chr[1],locus=i,Lead_snp=subset(dat1,P_LC==min(subset(dat1,cytoband==i)$P_LC))$vid,No_of_SNPs_tested=nrow(subset(dat1,cytoband==i)),No_of_Fvars=nrow(subset(dat1,cytoband==i & Type=="Fvars")),No_of_Fvars_with_annotations=nrow(subset(dat2,cytoband==i)),Fvars_with_annotations=paste(subset(dat2,cytoband==i)$vid,collapse=";"))
return(res)
}
i<- unique(dat2$cytoband)
cl <- makeCluster(10,type="FORK")
num <- parLapply(cl,i,func)
library(plyr)
results<-ldply(num,data.frame)
stopCluster(cl)
results=results[order(results$Chr),]

write.table(results,"fvars_region_anno_stat.txt",row.names=F,col.names=T,quote=F,sep="\t")

#fvars annos#
fvars_annos=dat2[,c("chr","pos","vid","ref","alt","region","cytoband","log2FoldChange_allele","AF_CN_LC","P_LC","A549_coreMarks","IMR90_coreMarks","NHLF_coreMarks","Fetallung_coreMarks","Lung_coreMarks","chromHMM_active_tag","EncodeLung_peak","EncodeLung_ATAC_Dnase_tag","EncodeLung_H3K27ac_H3K4me3_H3K4me1_H3K9ac_tag","NMU_Lung_eqtl_gene","TFBS_HOCOMOCO","TFBS_jASPAR","TFBS_SNP2TFBS")]
write.table(fvars_annos,"fvars_anno_stat.txt",row.names=F,col.names=T,quote=F,sep="\t")

for (i in unique(fvars_annos$cytoband)){
tmp=subset(fvars_annos,cytoband==i)[,c("chr","pos","cytoband","vid")]
write.table(tmp,paste(i,"_anno.fvars",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
}

#NMU_eGenes#
egenes=unique(data.frame(strsplit(paste0(fvars_annos$NMU_Lung_eqtl_gene,collapse=";"),";")))
names(egenes)="gene_id"
anno=fread("genecodev19_anno.txt")
anno=subset(anno,symbol%in%egenes$gene_id)[,"ENSG"]
names(anno)="gene_id"

write.table(anno,"NMU_Lung_eqtl_gene.ENSG",row.names=F,col.names=T,quote=F,sep="\t")
