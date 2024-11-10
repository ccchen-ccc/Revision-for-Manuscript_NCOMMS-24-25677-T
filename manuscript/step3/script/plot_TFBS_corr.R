setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step3/data")

library(data.table)
library(stringr)
library(ggcorrplot)
library(ggpubr)
library(stringr)

#SNP2TFBS#
SNP2TFBS=fread("results_SNP2TFBS.txt",header=F)
SNP2TFBS$vid=paste(substring(SNP2TFBS$V1,4,nchar(SNP2TFBS$V1)),SNP2TFBS$V3,SNP2TFBS$V4,SNP2TFBS$V5,sep=":")
SNP2TFBS$TF=sapply(strsplit(SNP2TFBS$V6,';'),'[',2)
SNP2TFBS$TF=substring(SNP2TFBS$TF,4,nchar(SNP2TFBS$TF))

SNP2TFBS$ScoreDiff=sapply(strsplit(SNP2TFBS$V6,';'),'[',3)
SNP2TFBS$ScoreDiff=substring(SNP2TFBS$ScoreDiff,11,nchar(SNP2TFBS$ScoreDiff))

SNP2TFBS$MATCH=sapply(strsplit(SNP2TFBS$V6,';'),'[',1)
SNP2TFBS$MATCH=as.numeric(substring(SNP2TFBS$MATCH,7,nchar(SNP2TFBS$MATCH)))

vid=unlist(sapply(1:length(SNP2TFBS$vid),function(x){rep(subset(SNP2TFBS,vid==SNP2TFBS$vid[x])$vid,subset(SNP2TFBS,vid==SNP2TFBS$vid[x])$MATCH)}))
TF=unlist(sapply(1:length(SNP2TFBS$vid),function(x){strsplit(subset(SNP2TFBS,vid==SNP2TFBS$vid[x])$TF,',')}))
ScoreDiff=as.numeric(unlist(sapply(1:length(SNP2TFBS$vid),function(x){strsplit(subset(SNP2TFBS,vid==SNP2TFBS$vid[x])$ScoreDiff,',')})))

SNP2TFBS_2=data.frame(vid=vid,TFBS_SNP2TFBS=TF,ScoreDiff=ScoreDiff)
SNP2TFBS_2$TFBS_SNP2TFBS=toupper(SNP2TFBS_2$TFBS_SNP2TFBS)
write.table(SNP2TFBS_2,"TFBS_SNP2TFBS.txt",sep="\t",quote=F,col.names=T,row.names=F)

SNP2TFBS_2$abs_ScoreDiff=abs(SNP2TFBS_2$ScoreDiff)
SNP2TFBS_2=SNP2TFBS_2[order(SNP2TFBS_2$vid,-SNP2TFBS_2$abs_ScoreDiff),]
SNP2TFBS_2=SNP2TFBS_2[!duplicated(SNP2TFBS_2$vid),]

mpra=fread("sig05_case_AD_SC_vid_alleleEffect_3cell_cov_rep6_0.01_filter.txt")
mpra=mpra[,c("vid","Fvars","log2FoldChange_allele","padj_allele")]

mpra=mpra[order(mpra$vid,-mpra$Fvars,mpra$padj_allele),]
mpra=mpra[!duplicated(mpra$vid),]
mpra$abs_log2FoldChange_allele=abs(mpra$log2FoldChange_allele)

dat2=merge(mpra,SNP2TFBS_2)

p1 = ggplot(subset(dat2,ScoreDiff!=0 & Fvars==1), aes(x = log2FoldChange_allele, y = ScoreDiff)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("MPRA log2FoldChange") + ylab("SNP2TFBS score")+
				stat_cor(method = 'pearson', aes(x = log2FoldChange_allele, y = ScoreDiff))

ggsave("corr_SNP2TFBS_sig-allele-log2FoldChange.pdf",p1,dpi=300)

p1 = ggplot(subset(dat2,ScoreDiff!=0), aes(x = log2FoldChange_allele, y = ScoreDiff)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("MPRA log2FoldChange") + ylab("SNP2TFBS score")+
				stat_cor(method = 'pearson', aes(x = log2FoldChange_allele, y = ScoreDiff))

ggsave("corr_SNP2TFBS_allele-log2FoldChange.pdf",p1,dpi=300)

