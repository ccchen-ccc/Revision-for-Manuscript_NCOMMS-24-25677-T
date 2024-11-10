library(data.table)
library(ggplot2)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)

setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step3/data")

dat=fread("sig05_case_AD_SC_vid_alleleEffect_3cell_cov_rep6_0.01_filter.txt")
dat1=dat[order(dat$vid,-dat$padj_allele_tag,-dat$padj_scram_tag,dat$padj_allele),]
dat1=dat1[!duplicated(dat1$vid),]

dat1$Type=ifelse(dat1$padj_allele_tag==1 & dat1$padj_scram_tag==1,"Fvars","Non-fvars")


#预测评分比较#
pred=unique(dat1[,c("vid","Type","log2FoldChange_allele","cadd","linsight","LungENN_v1_all_max")])
pred=as.data.frame(pred)
pred$cadd=abs(pred$cadd)
pred$linsight=abs(pred$linsight)
pred$LungENN_v1_all_max=abs(pred$LungENN_v1_all_max)
pred$log2FoldChange_allele=abs(pred$log2FoldChange_allele)

my_comparison=list(c("Fvars","Non-fvars"))
#cadd#
p_cadd=ggviolin(pred,x="Type",y="cadd",fill="Type",palette=c('npg'),xlab = FALSE,ylab = "CADD score",font.label=list(size=15,face="bold"),panel.labs=FALSE)+stat_compare_means(comparison=my_comparison)
ggsave(file="./fvars_LC/pred_score/violin_pred_score_fvars_cadd.pdf",p_cadd,dpi=300)
#corr_cadd#
p2_cadd = ggplot(data = pred, aes(x = log2FoldChange_allele, y = cadd)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("MPRA log2FoldChange") + ylab("cadd score")+
				stat_cor(method = 'pearson', aes(x = log2FoldChange_allele, y = cadd))

ggsave(file="./fvars_LC/pred_score/corrplot_pred_score_fvars_cadd.pdf",p2_cadd,dpi=300)

p3_cadd = ggplot(data = subset(pred,Type=="Fvars"), aes(x = log2FoldChange_allele, y = cadd)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("MPRA log2FoldChange") + ylab("cadd score")+
				stat_cor(method = 'pearson', aes(x = log2FoldChange_allele, y = cadd))
ggsave(file="./fvars_LC/pred_score/corrplot_pred_score_fvars_sig_cadd.pdf",p3_cadd,dpi=300)


#linsight#
p_linsight=ggviolin(pred,x="Type",y="linsight",fill="Type",palette=c('npg'),xlab = FALSE,ylab = "LINSIGHT score",font.label=list(size=15,face="bold"),panel.labs=FALSE)+stat_compare_means(comparison=my_comparison)
ggsave(file="./fvars_LC/pred_score/violin_pred_score_fvars_linsight.pdf",p_linsight,dpi=300)
#corr_linsight#
p2_linsight = ggplot(data = pred, aes(x = log2FoldChange_allele, y = linsight)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("MPRA log2FoldChange") + ylab("linsight score")+
				stat_cor(method = 'pearson', aes(x = log2FoldChange_allele, y = linsight))

ggsave(file="./fvars_LC/pred_score/corrplot_pred_score_fvars_linsight.pdf",p2_linsight,dpi=300)

p3_linsight = ggplot(data = subset(pred,Type=="Fvars"), aes(x = log2FoldChange_allele, y = linsight)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("MPRA log2FoldChange") + ylab("linsight score")+
				stat_cor(method = 'pearson', aes(x = log2FoldChange_allele, y = linsight))
ggsave(file="./fvars_LC/pred_score/corrplot_pred_score_fvars_sig_linsight.pdf",p3_linsight,dpi=300)



#violin_LungENN_v1_all_max#
p_lungENN=ggviolin(pred,x="Type",y="LungENN_v1_all_max",fill="Type",palette=c('npg'),xlab = FALSE,ylab = "lungENN score",font.label=list(size=20,face="bold"),panel.labs=FALSE)+stat_compare_means(comparison=my_comparison)
ggsave(file="./fvars_LC/pred_score/violin_pred_score_fvars_LungENN_v1_all_max.pdf",p_lungENN,dpi=300)

#corr_LungENN_v1_all_max#
p2_lungENN = ggplot(data = pred, aes(x = log2FoldChange_allele, y = LungENN_v1_all_max)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("MPRA log2FoldChange") + ylab("lungENN score")+
				stat_cor(method = 'pearson', aes(x = log2FoldChange_allele, y = LungENN_v1_all_max))

ggsave(file="./fvars_LC/pred_score/corrplot_pred_score_fvars_LungENN_v1_all_max.pdf",p2_lungENN,dpi=300)
p3_lungENN = ggplot(data = subset(pred,Type=="Fvars"), aes(x = log2FoldChange_allele, y = LungENN_v1_all_max)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("MPRA log2FoldChange") + ylab("lungENN score")+
				stat_cor(method = 'pearson', aes(x = log2FoldChange_allele, y = LungENN_v1_all_max))
ggsave(file="./fvars_LC/pred_score/corrplot_pred_score_fvars_sig_LungENN_v1_all_max.pdf",p3_lungENN,dpi=300)

