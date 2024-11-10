###fig2a
library(tidyverse)
library(ggplot2)
library(ggbreak)
setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step3/data")
df <- read.table('sig05_case_AD_SC_vid_alleleEffect_3cell_cov_rep6_0.01_filter.txt',header = T,sep = '\t')
dat2 <- df %>% 
  arrange(desc(vid), desc(padj_allele_tag), desc(padj_scram_tag), padj_allele) %>% 
  distinct(vid,.keep_all = T) %>% 
  select(vid,log2FoldChange_expr,padj_expr,padj_allele,log2FoldChange_allele)
cut_off_padj = 0.01
cut_off_logFC = 0
dat2$change_expr <-  ifelse(dat2$padj_expr < as.numeric(cut_off_padj) & dat2$log2FoldChange_expr >= cut_off_logFC, 'Active',
                            ifelse(dat2$padj_expr < as.numeric(cut_off_padj) & dat2$log2FoldChange_expr < cut_off_logFC,'Repressed','Stable'))
p <- ggplot(
  
  dat2, aes(x = log2FoldChange_expr, y = -log10(padj_expr), colour=change_expr)) +
  geom_point(alpha=0.8, size=3) +
  scale_color_manual(values=c("#ff7f24","#4d85bd", "#d2dae2"))+
  geom_hline(yintercept = -log10(as.numeric(cut_off_padj)),lty=2,col="grey60",lwd=0.5) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  xlim(-5,6.5)+
  ylim(0,200)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="top", 
        legend.title = element_blank())+
  scale_y_break(c(80,100),space = 0.2,scales = 0.4)
ggsave("fig2a.pdf",width = 7,height = 7,p,dpi=320,device = 'pdf')
#############fig2b
library(tidyverse)
library(patchwork)
exp_enrich <- read.table('enrich_exp_3cell_cov_rep6_0.01_filter.txt',sep = '\t',header = T)
exp_enrich <- exp_enrich[c(1,2,4:7,10,11),c('p_value_expr_pos_tag','OR_value_expr_pos_tag','p_value_expr_neg_tag','OR_value_expr_neg_tag')] %>% 
  mutate(name=c('DNase','ATAC','H3K27ac','H3K4me3','H3K4me1','H3K9ac','H3K36me3','H3K9me3'),
         OR_value_expr_pos_tag=str_replace_all(OR_value_expr_pos_tag,pattern = ' ',''),
         OR_value_expr_pos_tag=str_replace_all(OR_value_expr_pos_tag,pattern = '\\(','-'),
         OR_value_expr_pos_tag=str_replace_all(OR_value_expr_pos_tag,pattern = '\\)',''),
         OR_value_expr_neg_tag=str_replace_all(OR_value_expr_neg_tag,pattern = ' ',''),
         OR_value_expr_neg_tag=str_replace_all(OR_value_expr_neg_tag,pattern = '\\(','-'),
         OR_value_expr_neg_tag=str_replace_all(OR_value_expr_neg_tag,pattern = '\\)',''))
exp_enrich <- exp_enrich %>% 
  separate(OR_value_expr_pos_tag,into = c("OddsRatio_pos","OR_pos_low",'OR_pos_high'),sep = '-') %>% 
  separate(OR_value_expr_neg_tag,into = c("OddsRatio_neg","OR_neg_low",'OR_neg_high'),sep = '-')
pos <- exp_enrich[,c(1:4,9)]
pos <- pos %>% 
  mutate(across(.cols = 1:4, .fns = as.numeric),
         name=factor(name, levels = c('DNase','ATAC','H3K4me3','H3K9ac','H3K27ac','H3K4me1','H3K36me3','H3K9me3')))
neg <- exp_enrich[,c(5:9)]
neg <- neg %>% 
  mutate(across(.cols = 1:4, .fns = as.numeric),
         name=factor(name, levels = c('DNase','ATAC','H3K4me3','H3K9ac','H3K27ac','H3K4me1','H3K36me3','H3K9me3')))
pos <- pos %>% 
  mutate(sig=ifelse(p_value_expr_pos_tag<=0.001,'***',
                    ifelse(p_value_expr_pos_tag<=0.01&p_value_expr_pos_tag>=0.001,'**','')),
         group='Active')
neg <- neg %>% 
  mutate(sig=ifelse(p_value_expr_neg_tag<=0.001,'***',
                    ifelse(p_value_expr_neg_tag<=0.01&p_value_expr_neg_tag>=0.001,'**','')),
         group='Repressed')

p2<-ggplot(neg,aes(x= name,y =log2(OddsRatio_neg),fill = group))+geom_bar(stat = "identity",width = 0.7)+geom_errorbar(aes(ymin = log2(OR_neg_low), ymax = log2(OR_neg_high)),width = 0.2)+
  theme(axis.text.x = element_text(size = 10,angle = 60,vjust = 0.5),panel.grid = element_blank(),panel.background = element_rect(fill = "white"),axis.line =   element_line(color = "black",size = 0.5))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits = c(-2,1))+
  scale_fill_manual(values = "cyan4")+
  geom_text(aes(label = sig),position = position_dodge(width = 0.9), vjust = -2.5)+
  guides(fill= guide_legend(title = ""))
p2
p1<-ggplot(pos,aes(x= name,y = log2(OddsRatio_pos),fill = group))+geom_bar(stat = "identity",width = 0.7)+
  theme(axis.text.x = element_text(size = 10,angle = 90),panel.grid = element_blank(),panel.background = element_rect(fill = "white"),axis.line =   element_line(color = "black",size = 0.5))+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x  = element_blank()
  )+scale_fill_manual(values = "chocolate1")+
  geom_text(aes(label = sig),position = position_dodge(width = 0.9), vjust = -1)+
  scale_y_continuous(limits = c(-1,2))+geom_errorbar(aes(ymin = log2(OR_pos_low), ymax = log2(OR_pos_high)), width = 0.2)+
  guides(fill= guide_legend(title = "",labels = c("Active","Repressed"),values = c("chocolate1","cyan4")))
p1
p1/p2
ggsave(p1/p2,filename = 'fig2b.pdf',path = './',width = 8,height =6 ,device = 'pdf')
######fig2c
library(tidyverse)
library(ggplot2)
df <- read.table('fvars_region_stat.txt',sep = '\t',header = T)
df <- df %>% mutate(num=ifelse(No_of_Fvars>=5,5,No_of_Fvars))
df2 <- df %>% group_by(num) %>% 
  summarise(n=n())
df2 <- df2 %>% mutate(num=as.character(num)) %>% 
  mutate(num=ifelse(num=="5","5+",num))
ggplot(data=df2,aes(x=num,y=n))+
  geom_bar(stat = 'identity',fill='black',width=0.9)
ggsave(p1/p2,filename = 'fig2c.pdf',path = './',width = 8,height =6 )###保存后ai美化
##########fig2d
dat <- read.table('sig05_case_AD_SC_vid_alleleEffect_3cell_cov_rep6_0.01_filter.txt',header = T,sep = '\t')
library(tidyverse)
library(ggplot2)
dat2 <- dat %>% 
  arrange(desc(vid), desc(padj_allele_tag), desc(padj_scram_tag), padj_allele) %>% 
  distinct(vid,.keep_all = T) %>% 
  select(vid,log2FoldChange_expr,padj_expr,padj_allele,log2FoldChange_allele)
cut_off_padj = 0.01
cut_off_logFC = 0 
dat2$change_allele = ifelse(dat2$padj_allele < as.numeric(cut_off_padj) & dat2$log2FoldChange_allele >= cut_off_logFC, 'Active',
                            ifelse(dat2$padj_allele < as.numeric(cut_off_padj) & dat2$log2FoldChange_allele < cut_off_logFC,'Repressed','Stable'))
max(dat2$log2FoldChange_allele)
p <- ggplot(
  dat2, aes(x = log2FoldChange_allele, y = -log10(padj_allele), colour=change_allele)) +
  geom_point(alpha=0.8, size=3) +
  scale_color_manual(values=c("#ff7f24","#4d85bd", "#d2dae2"))+
  geom_hline(yintercept = -log10(as.numeric(cut_off_padj)),lty=2,col="grey60",lwd=0.5) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  xlim(-3,2)+
  ylim(0,300)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="top", 
        legend.title = element_blank())
ggsave("fig2d.pdf",device = "pdf",width = 7,height = 7,p,dpi=320)
##############################fig2e
library(UpSetR)
library(tidyverse)
anno <- read.table('sig05_case_AD_SC_vid_alleleEffect_3cell_cov_rep6_0.01_filter.txt',header = T,sep = '\t')
anno <- anno %>% 
  filter(Fvars==1) %>% 
  group_by(vid) %>% 
  slice_min(padj_allele) %>% 
  ungroup() %>% 
  mutate('Chromatin/Histone'=ifelse(EncodeLung_ATAC_Dnase_tag==1|EncodeLung_H3K27ac_H3K4me3_H3K4me1_H3K9ac_tag==1,vid,0),
         chromHMM=ifelse(chromHMM_active_tag==1,vid,0),
         TFBS=ifelse(is.na(TFBS_HOCOMOCO)&is.na(TFBS_jASPAR)&is.na(TFBS_SNP2TFBS),0,vid),
         eQTL=ifelse(is.na(NMU_Lung_eqtl_gene),0,vid)) %>% 
  select('Chromatin/Histone',chromHMM,TFBS,eQTL)
annote_list <- list(anno$'Chromatin/Histone'[anno$'Chromatin/Histone'!=0],
                    anno$chromHMM[anno$chromHMM!=0],
                    anno$TFBS[anno$TFBS!=0],
                    anno$eQTL[anno$eQTL!=0])
names(annote_list) <- c("Chromatin/Histone'","chromHMM","TFBS","eQTL")
#
pdf('./fig2e.pdf',width = 12,height = 10)
upset(fromList(annote_list),
      nsets = 200,
      nintersects = 90,
      order.by = "freq",
      keep.order = F,
      mb.ratio = c(0.6,0.4),
      text.scale = 5,
      point.size = 8,
      line.size = 4
)
dev.off()
##########fig2f
library(tidyverse)
library(ggplot2)
anno <- read.table('sig05_case_AD_SC_vid_alleleEffect_3cell_cov_rep6_0.01_filter.txt',header = T,sep = '\t')
anno <- anno %>% 
  filter(Fvars==1) %>% 
  group_by(vid) %>% 
  slice_min(padj_allele) %>% 
  ungroup() %>% 
  mutate(TFBS=ifelse(is.na(TFBS_HOCOMOCO)&is.na(TFBS_jASPAR),0,1),
         eQTL=ifelse(is.na(NMU_Lung_eqtl_gene),0,1),
         'Chromatin/Histone'=ifelse(EncodeLung_ATAC_Dnase_tag==1|EncodeLung_H3K27ac_H3K4me3_H3K4me1_H3K9ac_tag==1,1,0)) %>% 
  dplyr::select(chromHMM_active_tag,TFBS,eQTL,'Chromatin/Histone')
anno$feature <- rowSums(anno)
percent <- prop.table(table(anno$feature)) %>% as.numeric()
feature<- c(0:4) %>% as.character()
dat <- cbind(feature=feature,percent=as.numeric(round(percent*100,digits = 1))) %>% as.data.frame()
dat$percent <- as.numeric(dat$percent)
dat <- dat %>% 
  arrange(percent)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank()
  )
bp<- ggplot(dat, aes(x="", y=percent, fill=feature))+
  geom_bar(width = 1, stat = "identity",color='black',size=1.5)
pie <- bp + coord_polar("y", start=0)
p=pie + scale_fill_manual(values=color)+blank_theme+
  theme(axis.text.x=element_blank())+
  geom_text(aes(label = paste0(percent, "%")), position = position_stack(vjust=0.5),size=5)
ggsave(p,filename = 'fig2f.pdf',device = 'pdf',dpi = 320,width = 8,height = 8,path = './')
#########fig2g、2h
library(tidyverse)
library(ggplot2)
library(ggpubr)
dat <- read.table('sig05_case_AD_SC_vid_alleleEffect_3cell_cov_rep6_0.01_filter.txt',header = T,sep = '\t')
dat2=dat[order(dat$vid,-dat$padj_allele_tag,-dat$padj_scram_tag,dat$padj_allele),]
dat2=dat2[!duplicated(dat2$vid),]
dat2$Type=ifelse(dat2$padj_allele_tag==1 & dat2$padj_scram_tag==1,"Fvars","Non-fvars")
pred=unique(dat2[,c("vid","Type","log2FoldChange_allele","cadd","linsight","deephaem_all_max","deephaem_lung_max","LungENN_v1_all_max","LungENN_v2_all_max","phyloP46way")])
pred=as.data.frame(pred)
pred$cadd=abs(pred$cadd)
pred$linsight=abs(pred$linsight)
pred$deephaem_all_max=abs(pred$deephaem_all_max)
pred$deephaem_lung_max=abs(pred$deephaem_lung_max)
pred$LungENN_v1_all_max=abs(pred$LungENN_v1_all_max)
pred$LungENN_v2_all_max=abs(pred$LungENN_v2_all_max)
pred$log2FoldChange_allele=abs(pred$log2FoldChange_allele)
my_comparison=list(c("Non-fvars","Fvars"))
p_lungENN=ggviolin(pred,x="Type",y="LungENN_v1_all_max",fill="Type",
                   palette=c('#4DBBD5FF','#E64B35FF'),xlab = FALSE,ylab = "lungENN score",
                   font.label=list(size=20,face="bold"),
                   panel.labs=FALSE)+stat_compare_means(comparison=my_comparison)
ggsave(file="fig2g.pdf",width = 7,height = 7,device = 'pdf',p_lungENN,dpi=320)
p3_lungENN = ggplot(data = subset(pred,Type=="Fvars"), aes(x = log2FoldChange_allele, y = LungENN_v1_all_max)) + 
  geom_point() + geom_smooth(method = lm) + theme_classic2()+
  scale_color_manual(values = c('#FF7400')) + 
  xlab("MPRA log2FoldChange") + ylab("lungENN score")+
  stat_cor(method = 'pearson', aes(x = log2FoldChange_allele, y = LungENN_v1_all_max))
ggsave(file="fig2f.pdf",width = 7,height = 7,device = 'pdf',p3_lungENN,dpi=320)
######
