args <- commandArgs(TRUE)
setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/annotation/allele/fvars_LC/eqtl/coloc")

library(data.table)
library(ggplot2)
library(ggrepel)

args[1] <- "4q22.1"


eqtl=fread(paste0("./",args[1],"/",args[1],"_NMU_Lung_eqtl.allpairs.txt"))
names(eqtl)=c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se")
vid=fread(paste0("./",args[1],"/plink/",args[1],".vid"),header=F)
anno=fread("/public/home/ccc2018/data/database/ref/gencodev19/genecodev19_anno.txt")
anno=anno[,c("gene_id","symbol")]
eqtl=merge(eqtl,anno)

dat=fread("/public/home/ccc2018/data/project/20210304_WGS_phasing/minimac4/with_fam/results/output_cn/shapeit4_cn_with_fam/plink_asso/combined_6_rpWGS/covar1_source7_GC1.1/asso_combined_6_rpWGS_chr_all_case.firth.covar1_source7_GC1.1.adjusted.filter.output")
gwas_id=subset(dat,ID%in%vid$V1)
names(gwas_id)=c("variant_id","p_gwas")

for (gene in unique(eqtl$symbol)) {

gene <- "FAM13A"

dat2=subset(eqtl,symbol==gene & variant_id %in% gwas_id$variant_id )[,c("variant_id","pval_nominal")]

#gene_eqtl=subset(gene_eqtl,pval_nominal>5e-08)#6p22.2
names(dat2)=c("variant_id","p_eqtl")
dat3=merge(gwas_id,dat2, by="variant_id")

mpra=fread("/public/home/ccc2018/data/project/MPRA/LC/GWAS/annotation/allele/alleleEffect_LC/sig05_case_AD_SC_vid_alleleEffect_3cell_cov_rep6_0.01.txt")

fvars=unique(subset(mpra,Fvars==1)[,c("vid","cytoband","GC","report_loci")])

dat3$label=ifelse(dat3$variant_id %in% c("4:89885714:T:C" , "4:89860843:G:A"),2,ifelse(dat3$variant_id %in% mpra$vid,1,0))

write.table(dat3,"4q22.1_fam13A.asso.tsv",row.names=F,col.names=T,sep="\t",quote=F)
dat3=fread("4q22.1_fam13A.asso.tsv")

plot=dat3
plot$p_gwas=-log10(plot$p_gwas)
plot$p_eqtl=-log10(plot$p_eqtl)

plot$label=as.factor(plot$label)
plot<-data.frame(plot[order(plot[,c("label")]),])
plot$u<-1
plot$u<-as.factor(plot$u)

pval=cor.test(plot$p_gwas,plot$p_eqtl,method="pearson")$p.value
r=round(cor(plot$p_gwas,plot$p_eqtl,method="pearson"),2)


plot$rs_id <- ifelse( plot$label == 2 , plot$variant_id , NA )
plot$size_use <- ifelse( plot$label == 2 , 1.2 , ifelse(plot$label == 1 , 1.01 , 1) )
plot$alpha_use <- ifelse( plot$label == 2 , 1 , 0.9 )
color_use <- c('#376B6D', '#9F353A', '#BDC0BA')
names(color_use) <- c('1', '2', '0' )


pdf <- ggplot(plot) + 
  geom_point(aes(x = p_gwas, y = p_eqtl, color = label , size = size_use , alpha = alpha_use , label = rs_id)) +
  geom_text_repel(aes(x = p_gwas, y = p_eqtl, label = rs_id) , nudge_x = -0.4, direction = "y", hjust = "right" , size = 3 , data = plot) + 
  scale_color_manual( breaks = names(color_use) ,values = color_use ) +
  xlab("GWAS -log10 (P)") +
  ylab(paste(gene," eqtl -log10 (P)",sep="")) +
  theme_bw() +
  scale_size(limits=c(1,2)) +
  scale_alpha(limits=c(0,1)) +
  theme( 
  	panel.border = element_blank(),
	axis.line = element_line(colour = "black"),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
	legend.position = 'none',
	axis.title = element_text(size = 12, face = "bold", color = 'BLACK' ),
	axis.text = element_text(size = 12, color = 'BLACK'),
)

ggsave(file=paste0("./",args[1],"/","COLOC_",args[1],"_",gene,".pdf"),pdf,dpi=300,width=4,height=3)


}