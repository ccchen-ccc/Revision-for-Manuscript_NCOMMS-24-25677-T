setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step3/data")
library(tidyverse)
dat <- read.table('sig05_case_AD_SC_vid_alleleEffect_3cell_cov_rep6_0.01_filter.txt',header = T,sep = '\t')
colnames(dat)
###根据数据构建homer需要的bed格式数据
#pos 
oligo_all <-  data.frame(chr=paste0('chr',dat$chr),start=dat$pos,end=dat$pos,id=dat$uniqid,pos_tag=dat$padj_expr_pos_tag,neg_tag=dat$padj_expr_neg_tag,expr_tag=dat$padj_expr_tag)
oligo_expr <- oligo_all %>% 
  filter(pos_tag==1) %>% 
  select(chr,start,end) %>% 
  distinct() %>% 
  mutate(start=start-60,end=end+60)
oligo_expr$id <- paste(oligo_expr$chr,oligo_expr$start,oligo_expr$end,sep = '_')
oligo_expr$strand <- '.'
oligo_expr$score <- '.'
write.table(oligo_expr,'oligo_expr_pos.txt',col.names = F,sep = '\t',row.names = F,quote = F)
# neg
oligo_expr_neg <- oligo_all %>% 
  filter(neg_tag==1) %>% 
  select(chr,start,end) %>%
  distinct() %>% 
  mutate(start=start-60,end=end+60)
oligo_expr_neg$id <- paste(oligo_expr_neg$chr,oligo_expr_neg$start,oligo_expr_neg$end,sep = '_')
oligo_expr_neg$strand <- '.'
oligo_expr_neg$score <- '.'
write.table(oligo_expr_neg,'oligo_expr_neg.txt',col.names = F,sep = '\t',row.names = F,quote = F)
#background 剔除总体中非pos与neg的
oligo_expr_background <- oligo_all %>% 
  select(chr,start,end) %>% 
  distinct() %>% 
  mutate(start=start-60,end=end+60) %>% 
  mutate(id=paste(chr,start,end,sep='_'))
oligo_expr_background <- oligo_expr_background %>% 
  filter(id%in%a) %>% 
  mutate(strand='.',
         score='.')
write.table(oligo_expr_background,'oligo_expr_background.txt',col.names = F,sep = '\t',row.names = F,quote = F)